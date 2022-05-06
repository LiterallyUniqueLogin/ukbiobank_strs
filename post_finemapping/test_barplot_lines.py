#!/usr/bin/env python3

import os

import bokeh.plotting
import bokeh.io
import bokeh.models

ukb = os.environ['UKB']

cats = [('asdf', 'a'), ('asdf', 'b'), ('basdf', 'c'), ('basdf', 'd')]
fig = bokeh.plotting.figure(
    x_range = bokeh.models.FactorRange(*cats)
)
fig.vbar(x=cats, top=[1,2,1,2],width=0.9)
fig.line(x=[('asdf', 'a', .2), ('asdf', 'b')], y=[3,4], color='red')
fig.line(x=[('asdf', 'a'), ('basdf', 'c')], y=[3,3], color='red')

bokeh.io.export_png(fig, filename=f'{ukb}/test.png')

