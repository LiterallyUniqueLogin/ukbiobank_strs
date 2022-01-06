#!/usr/bin/env python3

import argparse
import gzip

import bokeh.io
import bokeh.models
import bokeh.plotting
import numpy as np
import polars as pl

def linear_int_interpolate(c1, c2, dist, lower_cutoff):
    c_new = []
    for coord1, coord2 in zip(c1, c2):
        c_new.append(coord1 + round((coord2 - coord1)*(dist-lower_cutoff)/(1-lower_cutoff)))
    return c_new

parser = argparse.ArgumentParser()
parser.add_argument('plink_ld_file')
parser.add_argument('bim_file')
parser.add_argument('outfname')
parser.add_argument('str_pos', type=int)
parser.add_argument('--savevars')
parser.add_argument('--loadvars')
args = parser.parse_args()

assert (args.savevars is not None) != (args.loadvars is not None)

with gzip.open(args.plink_ld_file, 'rt') as f:
    n_vars = len(next(f).split())

print('Loading lds ...', flush=True, end='')
with gzip.open(args.plink_ld_file) as f:
    lds = pl.read_csv(
        f,
        sep='\t',
        has_header=False,
        dtypes={f'column_{i}': float for i in range(1, n_vars+1)}
    ).to_numpy()
print(f' done. n_vars loaded: {lds.shape[0]}', flush=True)

poses = []
# otherwise assume that variants have already been filtered
if args.savevars is not None:
    closeness_cutoff = .4
    neighbors_cutoff = 20
    keep = []
    for i in range(n_vars):
        print(f"Working on keeps at row {i}", flush=True, end='\r')
        if (
            np.sum(lds[i, :i] >= closeness_cutoff) +
            np.sum(lds[(i+1):, i] >= closeness_cutoff)
            >= neighbors_cutoff
        ):
            keep.append(i)
    lds = lds[:, keep][keep, :]
    print('Done working on keeps. n_vars left={lds.shape[0]}', flush=True)

    with open(args.bim_file) as bim_file, \
            open(args.savevars, 'w') as out_variant_file:
        for line_num, line in enumerate(bim_file):
            split = line.split()
            if line_num in keep:
                poses.append(int(split[3]))
                out_variant_file.write('\t'.join((split[0], split[3], split[4], split[5])) + '\n')
else:
    with open(args.loadvars) as variant_file:
        variants = variant_file.readlines()
    with open(args.bim_file) as bim_file:
        keep = []
        for line_num, line in enumerate(bim_file):
            split = line.split()
            if ('\t'.join((split[0], split[3], split[4], split[5])) + '\n') in variants:
                poses.append(int(split[3]))
                keep.append(line_num)
    assert len(keep) == len(variants)
    lds = lds[:, keep][keep, :]
    print('Done loading variants.', flush=True)

assert len(poses) == len(keep)
poses = np.array(poses)
assert np.all((poses[1:] - poses[:-1]) >= 0)
#print(poses)

tick_step = 50000
next_tick_pos = ((poses[0]//tick_step)+1)*tick_step
minor_tick_locs = []
major_tick_locs = []
major_tick_labels = []
str_loc = None
for idx, pos in enumerate(poses):
    if pos > next_tick_pos:
        if next_tick_pos % (tick_step*2) == 0:
            major_tick_locs.append(idx-1)
            major_tick_labels.append(next_tick_pos)
        else:
            minor_tick_locs.append(idx-1)
        next_tick_pos += tick_step
    if pos > args.str_pos and str_loc is None:
        str_loc = idx - 1
assert str_loc is not None
#print(str_loc, minor_tick_locs, major_tick_locs, major_tick_labels)


'''
cutoff = .4
keep = []
for i in range(n_vars):
    print(f"Working on keeps at row {i}", flush=True, end='\r')
    if np.any(lds[i, :i] >= cutoff) or np.any(lds[(i+1):, i] >= cutoff):
        keep.append(i)
print('Done working on keeps.           ', flush=True)

lds = lds[:, keep][keep, :]
print(lds.shape)

set_join_cutoff = 0.7
set_count_cutoff = 40
sets = [{i} for i in range(n_vars)]
set_joins = np.where(lds > set_join_cutoff)
n_joins = len(set_joins[0])
for count, (i,j) in enumerate(zip(*set_joins)):
    print(f"Working on sets {count}/{n_joins}... ", end='\r', flush=True)
    sets[i].add(j)
    sets[j] = sets[i]
keep = set()
for set_ in sets:
    if len(set_) > set_count_cutoff:
        keep = keep.union(set_)
keep = list(keep)
lds = lds[:, keep][keep, :]
print('Done working on sets.                                    ', flush=True)
print(lds.shape)
'''

n_vars = lds.shape[0]

size = 1200
fig = bokeh.plotting.figure(
    title='LD heatmap',
    width=size,
    height=size//2
)

fig.background_fill_color = None
fig.border_fill_color = None
fig.grid.grid_line_color = None
fig.toolbar_location = None
fig.yaxis.visible = False
fig.xaxis.ticker = bokeh.models.FixedTicker(
    ticks=major_tick_locs,
    minor_ticks=minor_tick_locs
)
fig.xaxis.major_label_overrides = {
    loc: f'{pos:g}' for (loc, pos) in zip(major_tick_locs, major_tick_labels)
}

c1 = (235, 220, 52) #EBDC34
c2 = (235, 77, 42) #EB4D2A
widths = np.full(int(n_vars*(n_vars-1)/2), size/(n_vars-1))
heights = widths.copy()
xs = []
ys = []
colors = []
color_cutoff=.05
for i in range(n_vars):
    print(f'Working on values at row {i}', end='\r', flush=True)
    colors.extend([
        linear_int_interpolate(c1, c2, lds[i,j], color_cutoff)
        if lds[i, j] >= color_cutoff
        else (255, 255, 255)
        for j in range(i)
    ])
    xs.extend(np.arange(0.5*i, i, 0.5)*size/(n_vars-1))
    ys.extend(np.arange(i-1, -1, -1)*size/2/(n_vars-1) + size/2)# size/2+size/2*(i-1)/(n_vars-1), size/2-size/2/(n_vars-1), -size/2/(n_vars-1)))#  1/np.sqrt(2)*np.arange(size/2,size/2-i,-1))
print('Done working on values.           ', flush=True)
cds = bokeh.models.ColumnDataSource(dict(
    widths=widths,
    heights=heights,
    xs=xs,
    ys=ys,
    colors=colors
))
fig.rect(
    angle=45,
    angle_units='deg',
    x='xs',
    y='ys',
    height='heights',
    width='widths',
    color='colors',
    source = cds,
    dilate=True,
)
fig.star(
    x=[str_loc*size/(n_vars - 1)],
    y=[size/2],
    color='green',
    size=10
)
bokeh.io.export_png(fig, filename=args.outfname)


