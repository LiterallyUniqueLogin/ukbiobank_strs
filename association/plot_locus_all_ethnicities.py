#!/usr/bin/env python3

import argparse
import ast
import json
import os
import sys

import bokeh.io
import bokeh.layouts
import bokeh.models
import bokeh.plotting
import cyvcf2
import numpy as np
import polars as pl

import python_array_utils as utils

ukb = os.environ['UKB']
sys.path.insert(0, f'{ukb}/../trtools/repo')

import trtools.utils.tr_harmonizer as trh

def fix_cols(cols):
    col_sightings = {}
    for col in cols:
        if col not in col_sightings:
            col_sightings[col] = 1
            yield col
        else:
            yield col + '__' + str(col_sightings[col])
            col_sightings[col] += 1

parser = argparse.ArgumentParser()
parser.add_argument('outloc')
parser.add_argument('imputed_vcf')
parser.add_argument('assoc_results_json', help='json dict str enthnicity to assoc_result')
parser.add_argument('pheno_datas_json', help='json dict str enthnicity to pheno data')
parser.add_argument('chrom', type=int)
parser.add_argument('pos', type=int)
parser.add_argument('phenotype')
parser.add_argument('dosage_cutoff', type=float)
parser.add_argument('--unit')
parser.add_argument('--binary', action='store_true', default=False)
args = parser.parse_args()

assert bool(args.unit) or args.binary

if not args.binary:
    y_axis_label='Mean ' + args.phenotype.replace('_', ' ') + f' ({args.unit})'
else:
    y_axis_label='Fraction '+ args.phenotype.replace('_', ' ') + ' cases'

ethnicities = ['white_brits', 'black', 'south_asian', 'chinese', 'irish', 'white_other']
assoc_results_d = json.loads(args.assoc_results_json)
pheno_datas_d = json.loads(args.pheno_datas_json)
assert set(assoc_results_d.keys()) == set(pheno_datas_d.keys()) == set(ethnicities)

figs = []
for ethnicity in ethnicities:
    if not args.binary:
        stat_name = 'mean'
    else:
        stat_name = 'fraction'

    result = pl.scan_csv(
        assoc_results_d[ethnicity],
        sep='\t',
        dtypes={'locus_filtered': str}
    ).filter(
        (pl.col('chrom') == args.chrom) & (pl.col('pos') == args.pos)
    ).collect().select([ # have to collect first due to some sort of bug
        'motif',
        '0.05_significance_CI',
        '5e-8_significance_CI',
        f'{stat_name}_{args.phenotype}_per_single_dosage',
    ])
    assert result.shape[0] == 1
        
    pheno_data = np.load(pheno_datas_d[ethnicity])

    bgen_samples = []
    with open(f'{ukb}/microarray/ukb46122_hap_chr1_v2_s487314.sample') as samplefile:
        for num, line in enumerate(samplefile):
            if num <= 1:
                # skip first two lines
                continue
            bgen_samples.append(line.split()[0])
    assert len(bgen_samples) == 487409
    samples_array = np.array(bgen_samples, dtype=float).reshape(-1, 1)

    merged_arr = utils.merge_arrays(samples_array, pheno_data)
    unfiltered_subset = ~np.isnan(merged_arr[:, 1])
    n_samples = np.sum(unfiltered_subset)

    vcf = cyvcf2.VCF(args.imputed_vcf)
    found_rec = False
    for record in vcf(f'{args.chrom}:{args.pos}-{args.pos}'):
        if record.POS < args.pos:
            continue
        if record.INFO.get('PERIOD') is None:
            continue

        assert not found_rec
        found_rec = True

        trrecord = trh.HarmonizeRecord(vcfrecord=record, vcftype='beagle-hipstr')

        len_alleles = [trrecord.ref_allele_length] + trrecord.alt_allele_lengths
        len_alleles = [round(allele_len, 2) for allele_len in len_alleles]

        ap1 = trrecord.format['AP1']
        ap1 = np.concatenate((1 - np.sum(ap1, axis=1).reshape(-1, 1), ap1), axis=1)
        ap2 = trrecord.format['AP2']
        ap2 = np.concatenate((1 - np.sum(ap2, axis=1).reshape(-1, 1), ap2), axis=1)

    # TODO this needs better testing
    subset_summed_dosages = {}
    for aidx1, len_allele1 in enumerate(len_alleles):
        for aidx2, len_allele2 in enumerate(len_alleles):
            summed_len = len_allele1 + len_allele2
            if summed_len not in subset_summed_dosages:
                subset_summed_dosages[summed_len] = 0
            subset_summed_dosages[summed_len] += np.sum(np.multiply(
                ap1[unfiltered_subset, aidx1], ap2[unfiltered_subset, aidx2]
            ))

    assert np.isclose(sum(subset_summed_dosages.values()), n_samples)

    alleles = list(subset_summed_dosages.keys())
    alleles_copy = alleles.copy()
    for allele in alleles_copy:
        if subset_summed_dosages[allele] < args.dosage_cutoff:
            alleles.remove(allele)
    alleles = sorted(alleles)

    mean_per_dosage = {float(allele): val for allele, val in ast.literal_eval(result[f'{stat_name}_{args.phenotype}_per_single_dosage'].to_numpy()[0]).items()}
    ci5e_2 = {float(allele): val for allele, val in ast.literal_eval(result['0.05_significance_CI'].to_numpy()[0]).items()}
    ci5e_8 = {float(allele): val for allele, val in ast.literal_eval(result['5e-8_significance_CI'].to_numpy()[0]).items()}
    #y_min = min(ci5e_8[allele][0] for allele in alleles)
    #y_max = max(ci5e_8[allele][1] for allele in alleles)

    fig_kws = dict(
        width = 600,
        height = 600,
        y_axis_label=y_axis_label,
        x_axis_label='Sum of allele lengths (repeat copies)',
        title=f'{ethnicity} (n={n_samples})'
    )
    if figs:
        fig_kws['y_range'] = figs[0].y_range
    figure = bokeh.plotting.figure(**fig_kws)
    figs.append(figure)

    figure.grid.grid_line_color = None
    figure.background_fill_color = None
    figure.border_fill_color = None
    figure.toolbar_location = None
    figure.title.text_font_size = '18px'
    figure.axis.axis_label_text_font_size = '18px'
    figure.axis.major_label_text_font_size = '14px'

    figure.varea(alleles, [ci5e_2[allele][1] for allele in alleles], [ci5e_8[allele][1] for allele in alleles], color="red", alpha=0.2, legend_label='1 - 5e-8 Confidence Interval')
    figure.varea(alleles, [ci5e_2[allele][0] for allele in alleles], [ci5e_2[allele][1] for allele in alleles], color="red", alpha=0.4, legend_label='0.95 Confidence Interval')
    figure.varea(alleles, [ci5e_8[allele][0] for allele in alleles], [ci5e_2[allele][0] for allele in alleles], color="red", alpha=0.2)
    figure.line(alleles, [mean_per_dosage[allele] for allele in alleles], line_width=2, color="black")
    figure.circle(alleles, [mean_per_dosage[allele] for allele in alleles], color="black", size=6, legend_label='mean')
    figure.legend.label_text_font_size = '10px'

    #figure.y_range = bokeh.models.Range1d(y_min - 0.05*(y_max-y_min), y_max + 0.05*(y_max-y_min))


headers = []
headers.append(bokeh.models.Div(
    text=f'STR {args.chrom}:{args.pos} Repeat Unit:{result["motif"].to_numpy()[0]}',
    align="start",
    #text_font_size='18px',
    height=70,
    sizing_mode='stretch_width',
    style={'font-size': '400%'}
))
headers.append(bokeh.models.Div(
    text=args.phenotype.replace('_', ' ').capitalize() + " vs genotype",
    align="start",
    #text_font_size='18px',
    height=70,
    sizing_mode='stretch_width',
    style={'font-size': '400%'}
))

footers = []
footers.append(bokeh.models.Div(
    text="Phenotype values are unadjusted for covariates",
    align="end",
    height=20,
    sizing_mode='stretch_width',
    style={'font-size': '200%'}
))
footers.append(bokeh.models.Div(
    text="People contribute to each genotype based on their prob. of having that genotype",
    align="end",
    height=20,
    sizing_mode='stretch_width',
    style={'font-size': '200%'}
))
footers.append(bokeh.models.Div(
    text="Only considers tested individuals",
    align="end",
    height=20,
    sizing_mode='stretch_width',
    style={'font-size': '200%'}
))
footers.append(bokeh.models.Div(
    text=f"Genotypes with dosages less than {int(args.dosage_cutoff)} are omitted",
    align="end",
    height=20,
    sizing_mode='stretch_width',
    style={'font-size': '200%'}
))

total_fig = bokeh.layouts.column(*headers, bokeh.layouts.row(figs), *footers)

bokeh.io.export_svg(total_fig, filename=f'{args.outloc}.svg')
bokeh.io.export_png(total_fig, filename=f'{args.outloc}.png')
