#!/usr/bin/env python3

import os

import bokeh.io
import bokeh.layouts
import bokeh.models

import plot_locus

ukb = os.environ['UKB']

#CBL
for thresh in 0.01, 0.001:
    '''
    figs = []
    for subpop in 'CBL_hom_begin_C_T_snp', 'CBL_hom_not_begin_C_T_snp':
        for phenotype, unit in ('platelet_count', '10^9 cells/L'), ('platelet_crit', '%'):
            figs.append(plot_locus.generate_figure(
                f'{ukb}/association/subpop_spot_test/{subpop}/white_brits/{phenotype}/chr11_119077000.tab',
                f'{ukb}/traits/subpop_subset_transformed_phenotypes/{subpop}/white_brits/{phenotype}.npy',
                11,
                119077000,
                phenotype,
                thresh,
                unit,
                None,
                True
            ))

    figs[0].add_layout(bokeh.models.Title(text="Homozygous imperfect", align="center", text_font_size='20px'), 'left')
    figs[2].add_layout(bokeh.models.Title(text="Homozygous perfect", align="center", text_font_size='20px'), 'left')

    for i in range(1,4):
        figs[i].x_range = figs[0].x_range

    figs[2].y_range = figs[0].y_range
    figs[3].y_range = figs[1].y_range

    bokeh.io.export_png(bokeh.layouts.column([bokeh.layouts.row(figs[:2]), bokeh.layouts.row(figs[2:])]), filename=f'{ukb}/association/locus_plots/CBL_begin_C_T_snp_{thresh}.png')
    '''

    figs = []
    for subpop in 'slc2a2_hom_snp', 'slc2a2_at_least_one_snp', 'slc2a2_hom_not_snp':
        for phenotype, unit in [('total_bilirubin', 'umol/L')]:
            figs.append(plot_locus.generate_figure(
                f'{ukb}/association/subpop_spot_test/{subpop}/white_brits/{phenotype}/chr3_170727702.tab',
                f'{ukb}/traits/subpop_subset_transformed_phenotypes/{subpop}/white_brits/{phenotype}.npy',
                3,
                170727702,
                phenotype,
                thresh,
                unit,
                None,
                True
            ))

    figs[0].add_layout(bokeh.models.Title(text="Homozygous imperfect", align="center", text_font_size='20px'), 'left')
    figs[1].add_layout(bokeh.models.Title(text="Imperfect in one or both alleles", align="center", text_font_size='20px'), 'left')
    figs[2].add_layout(bokeh.models.Title(text="Homozygous perfect", align="center", text_font_size='20px'), 'left')

    for i in range(1,3):
        figs[i].x_range = figs[0].x_range
        figs[i].y_range = figs[0].y_range

    bokeh.io.export_png(bokeh.layouts.column(figs), filename=f'{ukb}/association/locus_plots/slc2a2_snp_{thresh}.png')


