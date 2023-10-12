#!/usr/bin/env python3

import argparse
import ast

import bokeh.io
import bokeh.models
import bokeh.palettes
import bokeh.plotting
import polars as pl
import scipy.interpolate

def mod1000(num, precision):
    num = float(('{:.' + str(precision) + 'g}').format(num))
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    return num, magnitude

# from https://stackoverflow.com/a/45846841
def human_format(num):
    # shouldn't be called if num < 1000
    # returns X.Ymag (always including Y, even if 0) if num is less than 10 over the nearest multiple of 1000
    # returns XXmag if num is between 10 and 100 over the nearest multiple of 1000
    # returns XXXmag if num is between 100 and 1000 over the nearest multiple of 1000
    mod_num, magnitude = mod1000(num, 2)
    if mod_num >= 100:
        mod_num, magnitude = mod1000(num, 3)

    print_num = '{:f}'.format(mod_num).rstrip('0').rstrip('.')
    if mod_num < 10 and '.' not in print_num:
        print_num = f'{print_num}.0'

    return '{}{}'.format(print_num, ['', 'k', 'm', 'b', 't'][magnitude])
    #return '{}{}'.format('{:f}'.format(num).rstrip('0').rstrip('.'), ['', 'k', 'm', 'b', 't'][magnitude])

def fix_header(header):
    def fix_header_helper(_):
        part1 = header.rpartition('0.05_significance_CI')
        fix1 = part1[0] + 'foo' + part1[2]
        part2 = fix1.rpartition('5e-8_significance_CI')
        fix2 = part2[0] + 'bar' + part2[2]
        return fix2.split('\t')
    return fix_header_helper

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('outloc')
    parser.add_argument('phenotype')
    parser.add_argument('--chrom')
    parser.add_argument('--pos', type=int)
    parser.add_argument('--dosage-threshold', type=int)
    parser.add_argument('--dosage-fraction-threshold', type=float)
    parser.add_argument('--assoc-results', nargs='*')
    parser.add_argument('--datas-per-length-sum', nargs='*')
    parser.add_argument('--group-names', nargs='*')
    parser.add_argument('--unit')
    parser.add_argument('--binary', action='store_true', default=False)
    parser.add_argument('--residual-phenos', default=False, action='store_true')

    # UKB script convention : total_subset_dosage_per_summed_gt
    # associaTR dosages: total_dosage_per_summed_length
    # associaTR hardcalls: sample_count_per_summed_length
    parser.add_argument('--total-column-name')
    args = parser.parse_args()
   
    assert bool(args.datas_per_length_sum) != bool(args.assoc_results)

    if bool(args.assoc_results):
        assert args.pos is not None and args.chrom is not None
        n_results = len(args.assoc_results)
    else:
        assert bool(args.datas_per_length_sum)
        n_results = len(args.datas_per_length_sum)

    assert 1 <= n_results
    assert bool(args.dosage_threshold) + bool(args.dosage_fraction_threshold) <= 1
    if n_results > 1:
        assert len(args.group_names) == n_results

    dosage_dicts = []
    mean_per_dosages = []
    ci5e_2s = []
    if bool(args.assoc_results):
        if not args.binary:
            stat_name = 'mean'
        else:
            stat_name = 'fraction'

        cols = pl.read_csv(
            args.assoc_results[0],
            sep='\t',
            dtypes={'locus_filtered': str},
            n_rows=1
        ).columns

        if not args.residual_phenos:
            # UKB script naming convention
            small_ci_col = 'summed_0.05_significance_CI'
            pheno_col = f'{stat_name}_{args.phenotype}_per_summed_gt'
            if small_ci_col not in cols:
                # associaTR convention
                small_ci_col = 'summed_length_0.05_alpha_CI'
                pheno_col = f'{stat_name}_{args.phenotype}_per_summed_length'
        else:
            # UKB script convention
            small_ci_col = 'res_per_sum_0.05_significance_CI'
            pheno_col = f'{stat_name}_{args.phenotype}_residual_per_summed_gt'
            if small_ci_col not in cols:
                # associaTR convention
                small_ci_col = 'summed_length_0.05_alpha_CI'
                pheno_col = f'{stat_name}_residual_{args.phenotype}_per_summed_length'

        for assoc_results_fname in args.assoc_results:
            result = pl.scan_csv(
                assoc_results_fname,
                sep='\t',
                dtypes={'locus_filtered': str},
            ).filter(
                (pl.col('chrom').cast(str) == args.chrom) & (pl.col('pos') == args.pos)
            ).collect().select([ # have to collect first due to some sort of bug
                small_ci_col,
                pheno_col,
                args.total_column_name
            ])
            assert result.shape[0] == 1

            dosage_dicts.append({
                float(allele): val for allele, val in
                ast.literal_eval(result[args.total_column_name].to_numpy()[0]).items()
            })
            mean_per_dosages.append({float(allele): val for allele, val in ast.literal_eval(result[pheno_col].to_numpy()[0]).items()})
            ci5e_2s.append({float(allele): val for allele, val in ast.literal_eval(result[small_ci_col].to_numpy()[0]).items()})
    else:
        for data_per_length_sum in args.datas_per_length_sum:
            df = pl.read_csv(data_per_length_sum, sep='\t').filter(pl.col('count') > 1)
            dosage_dicts.append(dict(zip(df['length_sum'], df['count'])))
            mean_per_dosages.append(dict(zip(df['length_sum'], df['mean'])))
            ci5e_2s.append({allele: ast.literal_eval(tup) for allele, tup in zip(df['length_sum'], df['conf_int'])})

    figure = generate_figure(
        dosage_dicts,
        mean_per_dosages,
        ci5e_2s,
        args.group_names,
        args.phenotype,
        args.dosage_fraction_threshold if args.dosage_fraction_threshold else args.dosage_threshold,
        bool(args.dosage_fraction_threshold),
        args.unit,
        args.binary,
        args.residual_phenos,
    )

    bokeh.io.export_svg(figure, filename=f'{args.outloc}.svg')
    bokeh.io.export_png(figure, filename=f'{args.outloc}.png')

def generate_figure(
    dosage_dicts,
    mean_per_dosages,
    ci5e_2s,
    group_names,
    phenotype,
    threshold,
    is_fraction_threshold,
    unit,
    binary,
    use_residual_phenos,
):
    if is_fraction_threshold:
        assert 0 <= threshold <= 1
    elif threshold is not None:
        assert threshold > 1
        assert isinstance(threshold, int)

    y_axis_label = phenotype.replace('_', ' ')
    if use_residual_phenos:
        y_axis_label = 'residual ' + y_axis_label
    if not binary and unit is not None:
        y_axis_label = y_axis_label + f' ({unit})'
    elif binary:
        y_axis_label = 'fraction '+ y_axis_label + ' cases'

    y_axis_label = y_axis_label[0].upper() + y_axis_label[1:]

    figure = bokeh.plotting.figure(
        width = 800,
        height = 900,
        y_axis_label=y_axis_label,
        x_axis_label='Sum of allele lengths (repeat copies)',
        output_backend = 'svg'
    )
    figure.grid.grid_line_color = None
    figure.background_fill_color = None
    figure.border_fill_color = None
    figure.toolbar_location = None
    figure.axis.axis_label_text_font_size = '36px'
    figure.axis.major_label_text_font_size = '30px'
    figure.axis.minor_tick_line_color = None

    if len(dosage_dicts) == 1:
        markers = ['circle']
    elif len(dosage_dicts) == 2:
        markers = ['circle', 'triangle']
    else:
        markers = ['circle']*len(dosage_dicts)
    y_mins = []
    y_maxs = []

    assert len(dosage_dicts) == len(mean_per_dosages) == len(ci5e_2s)
    for assoc_idx, (dosage_dict, mean_per_dosage, ci5e_2) in enumerate(zip(dosage_dicts, mean_per_dosages, ci5e_2s)):
        assert dosage_dict.keys() == mean_per_dosage.keys() == ci5e_2.keys()

        alleles = list(dosage_dict.keys())
        if threshold:
            denominator = 1 if not is_fraction_threshold else sum(dosage_dict.values())
            alleles_copy = alleles.copy()

            for allele in alleles_copy:
                if dosage_dict[allele]/denominator < threshold:
                    alleles.remove(allele)
        alleles = sorted(alleles)

        y_mins.append(min(ci5e_2[allele][0] for allele in alleles))
        y_maxs.append(max(ci5e_2[allele][1] for allele in alleles))

        if len(dosage_dicts) == 1:
            color = 'dimgrey'
        else:
            if len(dosage_dicts) < 5:
                color =  bokeh.palettes.Colorblind[4][assoc_idx] # if assoc_idx == 0 else 3]
            else:
                color = bokeh.palettes.Category20[len(dosage_dicts)][assoc_idx]
        small_ci_legend_label = '95% CI'
        if len(dosage_dicts) > 1:
            small_ci_legend_label = f'{group_names[assoc_idx]} ' + small_ci_legend_label

        #for y_vals, offset, color, bold in [(ci5e_2, .45, 'red', 2)]: #(ci5e_8, .3, 'orange', 1), (ci5e_2, .45, 'red', 2):
        y_vals = ci5e_2
        offset = 0.25
        alpha = .8 #0.5
        # vertical bar connecting the whiskers
        min_y = min(y_vals[allele][0] for allele in alleles)
        max_y = max(y_vals[allele][1] for allele in alleles)

        figure.multi_line(
            [[allele, allele] for allele in alleles],
            [[y_vals[allele][0], y_vals[allele][1]] for allele in alleles],
            color=color,
            line_width=1,
            legend_label=small_ci_legend_label,
            alpha=alpha
        )
        # whiskers
        for direction in 0, 1:
            figure.multi_line(
                [[allele-offset, allele+offset] for allele in alleles],
                [[y_vals[allele][direction], y_vals[allele][direction]] for allele in alleles],
                color=color,
                line_width=3,
            )

        # totals
        cds = bokeh.models.ColumnDataSource(dict(
            x=[allele - 0.025*(max(alleles) - min(alleles)) for allele in alleles],
            y=[y_vals[allele][1] + (max_y-min_y)*0.01 for allele in alleles],
            text=[(f'{dosage_dict[allele]:,.0f}' if dosage_dict[allele] < 1000 else human_format(dosage_dict[allele])) for allele in alleles]
        ))
        figure.add_layout(bokeh.models.LabelSet(x='x', y='y',
                                             #x_offset=.05,
                                             text='text', source=cds, text_font_size='18px', text_color='black'))

        scatter_label = 'mean'
        if len(dosage_dicts) > 1:
            scatter_label = f'{group_names[assoc_idx]}'
            figure.line(alleles, [mean_per_dosage[allele] for allele in alleles], line_width=2, color=color)
            #figure.line(alleles, [mean_per_dosage[allele] for allele in alleles], line_width=1, color='black')
        else:
            figure.line(alleles, [mean_per_dosage[allele] for allele in alleles], line_width=2, color='black')
            color='black'
        figure.scatter(alleles, [mean_per_dosage[allele] for allele in alleles], marker=markers[assoc_idx],
                       color=color,
                       #fill_color=color, line_color='black',
                       size=7, legend_label=scatter_label)
        '''
        if len(dosage_dicts) > 1:
            figure.legend[0].items.insert(len(figure.legend[0].items), figure.legend[0].items[-2])
            del figure.legend[0].items[-3]
        '''

    figure.legend.label_text_font_size = '30px'
    y_min = min(y_mins)
    y_max = max(y_maxs)
    figure.y_range = bokeh.models.Range1d(y_min - 0.05*(y_max-y_min), y_max + 0.05*(y_max-y_min))

    return figure

if __name__ == '__main__':
    main()
