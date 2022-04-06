#!/usr/bin/env python3

import os
import os.path

import bokeh.io
import bokeh.layouts
import bokeh.models
import bokeh.plotting
import numpy as np
import polars as pl
import statsmodels.api as sm

ukb = os.environ['UKB']


for out_exp, suffix, y_axis_name in [
    (lambda pheno: (-pl.col(f'p_{pheno}').log10()).alias('out'), '_p_val', '-log10 p-value'),
    (lambda pheno: (pl.col(f'coeff_{pheno}').abs()/pl.col(f'se_{pheno}')).pow(2).alias('out'), '_z2_score', '(z score)^2'),
]:
    figs = [[]]
    for pheno, chrom, pos in [
        ('urate', 1, 82910782),
        ('eosinophil_count', 1, 101632792),
        ('alanine_aminotransferase', 1, 16506924),
        ('triglycerides', 1, 63300856),
        ('shbg', 1, 107622328),
        ('white_blood_cell_count', 4, 74909874),
        ('calcium', 3, 122085027),
        ('cystatin_c', 4, 77388199),
        ('haemoglobin_concentration', 6, 26066710),
        ('red_blood_cell_count', 6, 135479912)
    ]:
        xs = []
        ys = []
        for size in 1, 2, 4, 8, 16, 32, 64, 128:
            size *= 1000
            for replicate in range(1, 11):
                fname = f'{ukb}/sample_qc/subpop_runs/white_brits_subset_{size}_{replicate}/white_brits/{pheno}/combined_unrelated.sample'
                # 47 are missing
                if not os.path.exists(fname):
                    continue
                with open(fname) as samps:
                    xs.append(len(samps.readlines()) - 1)
                ys.append(pl.read_csv(
                    f'{ukb}/association/subpop_spot_test/white_brits_subset_{size}_{replicate}/white_brits/{pheno}/chr{chrom}_{pos}.tab',
                    sep='\t'
                ).select(out_exp(pheno))['out'].to_numpy()[0])

        sm_xs = np.hstack((np.array(xs).reshape(-1, 1), np.ones((len(xs), 1))))
        params = sm.OLS(ys, sm_xs).fit().params
        slope = params[0]
        intercept = params[1]

        figure = bokeh.plotting.figure(
            width=1200,
            height=1200,
            title=f'{pheno} {chrom}:{pos}',
            y_axis_label=y_axis_name,
            x_axis_label='sample size',
        )
        figure.toolbar_location = None
        figure.grid.grid_line_color = None
        figure.background_fill_color = None
        figure.border_fill_color = None
        figure.circle(xs, ys)
        figure.line([0, max(xs)], [intercept, max(xs)*slope + intercept])
        if len(figs[-1]) == 3:
            figs.append([])
        figs[-1].append(figure)

    figs.insert(0, [bokeh.models.Div(
        text='p-val vs sample-size'
    )])

    bokeh.io.export_png(bokeh.layouts.layout(figs), filename=f'{ukb}/association/plots/subsample_replicates{suffix}.png')
    bokeh.io.export_svg(bokeh.layouts.layout(figs), filename=f'{ukb}/association/plots/subsample_replicates{suffix}.svg')

