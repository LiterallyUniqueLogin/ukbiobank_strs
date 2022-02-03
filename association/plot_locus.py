#!/usr/bin/env python3

import argparse
import ast
import os
import sys

import bokeh.io
import bokeh.models
import bokeh.plotting
import numpy as np
import polars as pl

import python_array_utils as utils
import sample_utils

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
parser.add_argument('assoc_results')
parser.add_argument('pheno_data')
parser.add_argument('chrom', type=int)
parser.add_argument('pos')
parser.add_argument('phenotype')
parser.add_argument('dosage_fraction_threshold', type=float)
parser.add_argument('--unit')
parser.add_argument('--binary', action='store_true', default=False)
args = parser.parse_args()

assert bool(args.unit) or args.binary

assert 0 <= args.dosage_fraction_threshold <= 1

if not args.binary:
    y_axis_label='Mean ' + args.phenotype.replace('_', ' ') + f' ({args.unit})'
else:
    y_axis_label='Fraction '+ args.phenotype.replace('_', ' ') + ' cases'

figure = bokeh.plotting.figure(
    width = 600,
    height = 600,
    y_axis_label=y_axis_label,
    x_axis_label='Sum of allele lengths (repeat copies)'
)
figure.grid.grid_line_color = None
figure.background_fill_color = None
figure.border_fill_color = None
figure.toolbar_location = None
figure.title.text_font_size = '18px'
figure.axis.axis_label_text_font_size = '18px'
figure.axis.major_label_text_font_size = '14px'

if not args.binary:
    stat_name = 'mean'
else:
    stat_name = 'fraction'

result = pl.scan_csv(
    args.assoc_results,
    sep='\t',
    dtypes={'locus_filtered': str}
).filter(
    (pl.col('chrom') == args.chrom) & (pl.col('pos') == args.pos)
).collect().select([ # have to collect first due to some sort of bug
    'motif',
    '0.05_significance_CI',
    '5e-8_significance_CI',
    f'{stat_name}_{args.phenotype}_per_single_dosage',
    'total_subset_dosage_per_summed_gt'
])
assert result.shape[0] == 1
    
pheno_data = np.load(args.pheno_data)

bgen_samples = sample_utils.get_all_samples()
assert len(bgen_samples) == 487409
samples_array = np.array(bgen_samples, dtype=float).reshape(-1, 1)

merged_arr = utils.merge_arrays(samples_array, pheno_data)
unfiltered_subset = ~np.isnan(merged_arr[:, 1])
n_samples = np.sum(unfiltered_subset)

subset_summed_dosage_fractions = {
    float(allele): val for allele, val in
    ast.literal_eval(result['total_subset_dosage_per_summed_gt'].to_numpy()[0]).items()
}
total_dosage = np.sum(list(subset_summed_dosage_fractions.values()))
subset_summed_dosage_fractions = {key: val/total_dosage for key, val in subset_summed_dosage_fractions.items()}

alleles = list(subset_summed_dosage_fractions.keys())
alleles_copy = alleles.copy()
for allele in alleles_copy:
    if subset_summed_dosage_fractions[allele] < args.dosage_fraction_threshold:
        alleles.remove(allele)
alleles = sorted(alleles)

mean_per_dosage = {float(allele): val for allele, val in ast.literal_eval(result[f'{stat_name}_{args.phenotype}_per_single_dosage'].to_numpy()[0]).items()}
ci5e_2 = {float(allele): val for allele, val in ast.literal_eval(result['0.05_significance_CI'].to_numpy()[0]).items()}
ci5e_8 = {float(allele): val for allele, val in ast.literal_eval(result['5e-8_significance_CI'].to_numpy()[0]).items()}
y_min = min(ci5e_8[allele][0] for allele in alleles)
y_max = max(ci5e_8[allele][1] for allele in alleles)

figure.varea(alleles, [ci5e_2[allele][1] for allele in alleles], [ci5e_8[allele][1] for allele in alleles], color="red", alpha=0.2, legend_label='1 - 5e-8 Confidence Interval')
figure.varea(alleles, [ci5e_2[allele][0] for allele in alleles], [ci5e_2[allele][1] for allele in alleles], color="red", alpha=0.4, legend_label='0.95 Confidence Interval')
figure.varea(alleles, [ci5e_8[allele][0] for allele in alleles], [ci5e_2[allele][0] for allele in alleles], color="red", alpha=0.2)
figure.line(alleles, [mean_per_dosage[allele] for allele in alleles], line_width=2, color="black")
figure.circle(alleles, [mean_per_dosage[allele] for allele in alleles], color="black", size=6, legend_label='mean')
figure.legend.label_text_font_size = '10px'

figure.y_range = bokeh.models.Range1d(y_min - 0.05*(y_max-y_min), y_max + 0.05*(y_max-y_min))

figure.add_layout(
    bokeh.models.Title(text=f'STR {args.chrom}:{args.pos} Repeat Unit:{result["motif"].to_numpy()[0]}', align="center", text_font_size='18px'), "above"
)
figure.add_layout(
    bokeh.models.Title(text=args.phenotype.replace('_', ' ').capitalize() + " vs genotype", align="center", text_font_size='18px'), "above"
)

figure.add_layout(bokeh.models.Title(text="Phenotype values are unadjusted for covariates", align="center"), "below")
figure.add_layout(bokeh.models.Title(text="People contribute to each genotype based on their prob. of having that genotype", align="center"), "below")
figure.add_layout(bokeh.models.Title(text="Only considers tested individuals", align="center"), "below")
figure.add_layout(bokeh.models.Title(text=f"Genotypes with dosages less than {100*args.dosage_fraction_threshold}% of the population are omitted", align="center"), "below")

bokeh.io.export_svg(figure, filename=f'{args.outloc}.svg')
bokeh.io.export_png(figure, filename=f'{args.outloc}.png')
