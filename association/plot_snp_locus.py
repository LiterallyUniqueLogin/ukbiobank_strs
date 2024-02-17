import argparse

import bokeh.io
import bokeh.models
import bokeh.plotting
import bgen_reader
import numpy as np
import polars as pl
import statsmodels.stats.weightstats

import weighted_binom_conf

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

    print_num = '{:f}'.format(mod_num).rstrip('0').rstrip('.') # W: Formatting a regular string which could be an f-string
    if mod_num < 10 and '.' not in print_num:
        print_num = f'{print_num}.0'

    return '{}{}'.format(print_num, ['', 'k', 'm', 'b', 't'][magnitude]) # W: Formatting a regular string which could be an f-string
    #return '{}{}'.format('{:f}'.format(num).rstrip('0').rstrip('.'), ['', 'k', 'm', 'b', 't'][magnitude])

parser = argparse.ArgumentParser()
parser.add_argument('outprefix')
parser.add_argument('bgen')
parser.add_argument('pos', type=int)
parser.add_argument('ref')
parser.add_argument('alt')
parser.add_argument('bgen_samples') # $UKB/array_imputed/ukb46122_imp_chr1_v3_s487283.sample
parser.add_argument('pheno_file') # first column must be sample number
parser.add_argument('pheno_column')
parser.add_argument('y_label')
parser.add_argument('--binary', default=False, action='store_true')
args = parser.parse_args()

reader = bgen_reader.open_bgen(args.bgen, verbose=False)

idx = (reader.positions == args.pos) & (reader.allele_ids == f'{args.ref},{args.alt}')
assert sum(idx) == 1
idx = np.where(idx)[0][0]

gts = reader.read(idx).squeeze()
bgen_samples = pl.read_csv(args.bgen_samples, sep=' ').filter(pl.col('ID_1') != 0).rename({'ID_1': 'id'}).select('id')
gts_df = pl.DataFrame({'id': bgen_samples['id'], 'p_hom_ref': gts[:, 0], 'p_het': gts[:, 1], 'p_hom_alt': gts[:, 2]})

pheno_df = pl.read_csv(args.pheno_file, sep='\t')
pheno_df = pheno_df.rename({pheno_df.columns[0]: 'id', args.pheno_column: 'pheno'}).filter(~pl.col('pheno').is_null()).select(['id', 'pheno'])

df = gts_df.join(pheno_df, how='inner', on='id')
n_samples = df.shape[0]

means = []
ci_lowers = []
ci_uppers = []
for col in 'p_hom_ref', 'p_het', 'p_hom_alt':
    if df.filter(pl.col(col) != 0).select('pheno').unique().shape[0] < 1:
        means.append(np.nan)
        ci_lowers.append(np.nan)
        ci_uppers.append(np.nan)
        continue
    if df.filter(pl.col(col) != 0).select('pheno').unique().shape[0] == 1:
        means.append(df.filter(pl.col(col) != 0).select('pheno').unique()[0, 0])
        ci_lowers.append(np.nan)
        ci_uppers.append(np.nan)
        continue
    gts = df[col].to_numpy().flatten()
    phenos = df['pheno'].to_numpy().flatten()
    if not args.binary:
        mean_stats = statsmodels.stats.weightstats.DescrStatsW(phenos, weights=gts)
        mean = mean_stats.mean
        ci_lower, ci_upper = mean_stats.tconfint_mean()
    else:
        mean, ci_lower, ci_upper = weighted_binom_conf.weighted_binom_conf(gts, phenos, 0.95)
    means.append(mean)
    ci_lowers.append(ci_lower)
    ci_uppers.append(ci_upper)
ci_xs = [x for x in range(3) if ~np.isnan(ci_lowers[x])]
mean_xs = [x for x in range(3) if ~np.isnan(means[x])]

figure = bokeh.plotting.figure(
	width = 800,
	height = 900,
	y_axis_label=args.y_label,
	x_axis_label='# alt alleles',
	output_backend = 'svg'
)
figure.grid.grid_line_color = None
figure.background_fill_color = None
figure.border_fill_color = None
figure.toolbar_location = None
figure.axis.axis_label_text_font_size = '36px'
figure.axis.major_label_text_font_size = '30px'
figure.axis.minor_tick_line_color = None

figure.xaxis.ticker = [0, 1, 2]
figure.x_range = bokeh.models.Range1d(-0.5, 2.5)

y_max = max(ci_uppers[x] for x in ci_xs)
y_min = min(ci_lowers[x] for x in ci_xs)
figure.y_range = bokeh.models.Range1d(y_min - 0.05*(y_max-y_min), y_max + 0.05*(y_max-y_min))

color = 'dimgrey'

# vertical bar connecting the whiskers
alpha = .8
figure.multi_line(
    [(x, x) for x in ci_xs],
    [[ci_lowers[x], ci_uppers[x]] for x in ci_xs],
    color=color,
    line_width=1,
    legend_label='95% CI',
    alpha=alpha
)

# whiskers
offset = 0.25
for cis in ci_lowers, ci_uppers:
    figure.multi_line(
        [(x - offset, x + offset) for x in ci_xs],
        [[cis[x], cis[x]] for x in ci_xs],
        color=color,
        line_width=3,
    )

# total labels
cds = bokeh.models.ColumnDataSource(dict(
    x=[x - 0.05 for x in range(3)],
    y=[(ci_uppers[x] + (y_max-y_min)*0.01 if ~np.isnan(means[x]) else ci_uppers[1] + (y_max-y_min)*0.01) for x in range(3)],
    text=[f'{human_format(df.select(pl.col(col).sum())[0, 0])}, {df.select(pl.col(col).sum())[0, 0]*100/n_samples:.2g}%' for col in ['p_hom_ref', 'p_het', 'p_hom_alt']]
))
figure.add_layout(bokeh.models.LabelSet(
    x='x', y='y',
    #x_offset=.05,
    text='text', source=cds, text_font_size='18px', text_color='black'
))

figure.line(mean_xs, [means[x] for x in mean_xs], line_width=2, color='black')
figure.scatter(
    mean_xs,
    [means[x] for x in mean_xs],
    marker='circle',
    color='black',
    size=7, legend_label='mean'
)

figure.legend.label_text_font_size = '30px'
bokeh.io.export_png(figure, filename=f'{args.outprefix}.png')

