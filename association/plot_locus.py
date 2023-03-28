#!/usr/bin/env python3

import argparse
import ast

import bokeh.io
import bokeh.models
import bokeh.palettes
import bokeh.plotting
import numpy as np
import polars as pl
import scipy.interpolate

# from https://stackoverflow.com/a/45846841
def human_format(num):
    num = float('{:.2g}'.format(num))
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    return '{}{}'.format('{:f}'.format(num).rstrip('0').rstrip('.'), ['', 'k', 'm', 'b', 't'][magnitude])

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
    parser.add_argument('chrom', type=int)
    parser.add_argument('pos', type=int)
    parser.add_argument('phenotype')
    parser.add_argument('--dosage-threshold', type=int)
    parser.add_argument('--dosage-fraction-threshold', type=float)
    parser.add_argument('--assoc-results', nargs='+')
    parser.add_argument('--group-names', nargs='+', default=[])
    parser.add_argument('--unit')
    parser.add_argument('--binary', action='store_true', default=False)
    parser.add_argument('--residual-phenos', default=False, action='store_true')
    parser.add_argument('--publication', default=False, action='store_true')
    parser.add_argument('--count-per-length-sum', nargs='+')
    parser.add_argument('--stat-per-length-sum', nargs='+')
    parser.add_argument('--05-CI-per-length-sum', nargs='+')
    args = parser.parse_args()
   
    assert bool(args.total_subset_dosage_per_summed_gt) != bool(args.assoc_results)
    assert bool(args.total_subset_dosage_per_summed_gt) == bool(args.stat_per_summed_genotype) == bool(args.summed_05_significance_CI)

    if bool(args.assoc_results):
        n_results = len(args.assoc_results)
    else:
        assert bool(args.total_subset_dosage_per_summed_gt)
        n_results = len(args.total_subset_dosage_per_summed_gt)

    assert 1 >= n_results
    assert bool(args.dosage_threshold) + bool(args.dosage_fraction_threshold) == 1
    if n_results > 1:
        assert len(args.group_names) == len(args.assoc_results)
        assert args.dosage_threshold is not None

    if bool(args.assoc_results):
        dosage_test_dicts = []
        mean_per_dosages = []
        ci5e_2s = []
        if not args.binary:
            stat_name = 'mean'
        else:
            stat_name = 'fraction'

        if not args.residual_phenos:
            small_ci_col = 'summed_0.05_significance_CI'
            pheno_col = f'{stat_name}_{args.phenotype}_per_summed_genotype'
        else:
            small_ci_col = 'res_per_sum_0.05_significance_CI'
            pheno_col = f'{stat_name}_{args.phenotype}_residual_per_summed_gt'

        for assoc_results_fname in args.assoc_results:
            with open(assoc_results_fname) as tsv:
                header = tsv.readline().strip()
            result = pl.scan_csv(
                assoc_results_fname,
                sep='\t',
                dtypes={'locus_filtered': str},
                skip_rows=1,
                has_header=False,
                with_column_names = fix_header(header)
            ).filter(
                (pl.col('chrom') == args.chrom) & (pl.col('pos') == args.pos)
            ).collect().select([ # have to collect first due to some sort of bug
                small_ci_col,
                pheno_col,
                'total_subset_dosage_per_summed_gt'
            ])
            assert result.shape[0] == 1
           
            dosage_dict = {
                float(allele): val for allele, val in
                ast.literal_eval(result['total_subset_dosage_per_summed_gt'].to_numpy()[0]).items()
            }
            # TODO need to pass non fraction through regardless
            if bool(args.dosage_fraction_threshold):
                total_dosage = np.sum(list(dosage_dict.values()))
                dosage_test_dicts.append({key: val/total_dosage for key, val in dosage_dict.items()})
            else:
                dosage_test_dicts.append(dosage_dict)

            mean_per_dosages.append({float(allele): val for allele, val in ast.literal_eval(result[pheno_col].to_numpy()[0]).items()})
            ci5e_2s.append({float(allele): val for allele, val in ast.literal_eval(result[small_ci_col].to_numpy()[0]).items()})
    else:
        pass
#        parser.add_argument('--count-per-length-sum', nargs='+')
#        parser.add_argument('--stat-per-length-sum', nargs='+')
#        parser.add_argument('--05-CI-per-length-sum', nargs='+')
#        dosage_test_dicts = args.total_subset_dosage_per_summed_gt
#        mean_per_dosages = args.stat_per_summed_genotype
#        ci5e_2s = args.summed_05_significance_CI
     
    figure = generate_figure(
        dosage_test_dicts,
        mean_per_dosages,
        ci5e_2s,
        args.group_names,
        args.chrom,
        args.pos,
        args.phenotype,
        args.dosage_fraction_threshold if args.dosage_fraction_threshold else args.dosage_threshold,
        bool(args.dosage_fraction_threshold),
        args.unit,
        args.binary,
        args.residual_phenos,
        args.publication,
    )

    bokeh.io.export_svg(figure, filename=f'{args.outloc}.svg')
    bokeh.io.export_png(figure, filename=f'{args.outloc}.png')

def generate_figure(
    dosage_test_dicts,
    mean_per_dosages,
    ci5e_2s,
    group_names,
    chrom,
    pos,
    phenotype,
    threshold,
    is_fraction_threshold,
    unit,
    binary,
    use_residual_phenos,
    publication,
):
    assert bool(unit) or binary

    if is_fraction_threshold:
        assert 0 <= threshold <= 1
    else:
        assert threshold > 1
        assert isinstance(threshold, int)

    y_axis_label = phenotype.replace('_', ' ')
    if use_residual_phenos:
        y_axis_label = 'residual ' + y_axis_label
    if not binary:
        y_axis_label='mean ' + y_axis_label + f' ({unit})'
    else:
        # TODO check this
        y_axis_label='fraction '+ y_axis_label + ' cases'

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

    #markers = ['circle', 'triangle']
    markers = ['circle']*len(dosage_test_dicts)
    y_mins = []
    y_maxs = []

    assert len(dosage_test_dicts) == len(mean_per_dosages) == len(ci5e_2s)
    for assoc_idx, (dosage_test_dict, mean_per_dosage, ci5e_2) in enumerate(zip(dosage_test_dicts, mean_per_dosage, ci5e_2)):
        assert dosage_test_dict.keys() == mean_per_dosage.keys() == ci5e_2.keys()

        alleles = list(dosage_test_dict.keys())
        alleles_copy = alleles.copy()

        for allele in alleles_copy:
            if dosage_test_dict[allele] < threshold:
                alleles.remove(allele)
        alleles = sorted(alleles)

        y_mins.append(min(ci5e_2[allele][0] for allele in alleles))
        y_maxs.append(max(ci5e_2[allele][1] for allele in alleles))

        if len(dosage_test_dicts) == 1:
            color = 'dimgrey'
        else:
            if len(dosage_test_dicts) < 5:
                color =  bokeh.palettes.Colorblind[4][assoc_idx] # if assoc_idx == 0 else 3]
            else:
                color = bokeh.palettes.Category20[len(dosage_test_dicts)][assoc_idx]
        small_ci_legend_label = '95% CI'
        if len(dosage_test_dicts) > 1:
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
        if len(dosage_test_dicts) > 1:
            scatter_label = f'{group_names[assoc_idx]}'
            figure.line(alleles, [mean_per_dosage[allele] for allele in alleles], line_width=2, color=color)
            figure.line(alleles, [mean_per_dosage[allele] for allele in alleles], line_width=1, color='black')
        else:
            figure.line(alleles, [mean_per_dosage[allele] for allele in alleles], line_width=2, color='black')
            color='black'
        figure.scatter(alleles, [mean_per_dosage[allele] for allele in alleles], marker=markers[assoc_idx],
                       color=color,
                       #fill_color=color, line_color='black',
                       size=7, legend_label=scatter_label)
        '''
        if len(dosage_test_dicts) > 1:
            figure.legend[0].items.insert(len(figure.legend[0].items), figure.legend[0].items[-2])
            del figure.legend[0].items[-3]
        '''

    figure.legend.label_text_font_size = '30px'
    y_min = min(y_mins)
    y_max = max(y_maxs)
    figure.y_range = bokeh.models.Range1d(y_min - 0.05*(y_max-y_min), y_max + 0.05*(y_max-y_min))

    if not publication:
        figure.add_layout(
            bokeh.models.Title(text=f'STR {chrom}:{pos}', align="center", text_font_size='18px'), "above"
        )
        figure.add_layout(
            bokeh.models.Title(text=phenotype.replace('_', ' ').capitalize() + " vs genotype", align="center", text_font_size='18px'), "above"
        )

        figure.add_layout(bokeh.models.Title(text="Phenotype values are unadjusted for covariates", align="center"), "below")
        figure.add_layout(bokeh.models.Title(text="People contribute to each genotype based on their prob. of having that genotype", align="center"), "below")
        figure.add_layout(bokeh.models.Title(text="Only considers tested individuals", align="center"), "below")
        figure.add_layout(bokeh.models.Title(text=f"Genotypes with dosages less than {100*dosage_fraction_threshold}% of the population are omitted", align="center"), "below")

    return figure

if __name__ == '__main__':
    main()
