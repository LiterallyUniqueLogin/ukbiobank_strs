#!/usr/bin/env python3

import argparse

import bokeh.io
import bokeh.models
import bokeh.palettes
import bokeh.plotting
import numpy as np
import polars as pl

parser = argparse.ArgumentParser()
parser.add_argument('effect_sizes', nargs=4)
args = parser.parse_args()

fig = bokeh.plotting.figure(
    x_axis_label='effect size',
    y_axis_label='cumulative probability',
    width=1200,
    height=800,
    toolbar_location=None,
    output_backend='svg'
)

fig.grid.grid_line_color = None
fig.axis.axis_label_text_font_size = '36px'
fig.axis.major_label_text_font_size = '30px'
fig.title.text_font_size = '36px'

colors = bokeh.palettes.Colorblind4
legends = ('0.01% - 0.1%', '0.1% - 1%', '1% - 10%', '10% - 50%')
effect_sizess = []
for effect_size_fname in args.effect_sizes:
    with open(effect_size_fname) as effect_size_file:
        effect_sizess.append(sorted(abs(float(line)) for line in effect_size_file.readlines()))
max_effect_size = max([max(effect_sizes) for effect_sizes in effect_sizess])
for color, legend, effect_sizes in zip(colors, legends, effect_sizess):
    n = len(effect_sizes)
    effect_sizes = [0] + effect_sizes + [max_effect_size + 0.02]
    fig.step(
        effect_sizes,
        [0] + list(np.arange(n)/n + 1/n) + [1],
        mode = "after",
        legend_label = legend,
        color = color,
        line_width = 3
    )
fig.legend.location = 'center_right'
fig.legend.label_text_font_size = '30px'

bokeh.io.export_svg(fig, filename='bin_effect_sizes.svg')
bokeh.io.export_png(fig, filename='bin_effect_sizes.png')
