#!/usr/bin/env python3

import argparse

import bokeh.plotting
import bokeh.models
import polars as pl

import graphing_utils

parser = argparse.ArgumentParser()
parser.add_argument('out')
parser.add_argument('causal_candidates_table')
args = parser.parse_args()

candidates = pl.read_csv(args.causal_candidates_table, sep='\t').distinct(subset=['chrom', 'start_pos'])

n_intergenic = candidates.select((pl.col('relation_to_gene') == 'intergenic').sum().alias('out'))['out'].to_numpy()[0]
assert n_intergenic == candidates.select(pl.col('relation_to_gene').str.contains('intergenic').sum().alias('out'))['out'].to_numpy()[0]
n_nonconding = candidates.select(
    (
        (pl.col('relation_to_gene') != 'intergenic') &
        ~pl.col('relation_to_gene').str.contains('protein_coding')
    ).sum().alias('out')
)['out'].to_numpy()[0]

n_noncoding = candidates.select(
    (
        (pl.col('relation_to_gene') != 'intergenic') &
        ~pl.col('relation_to_gene').str.contains('protein_coding')
    ).sum().alias('out')
)['out'].to_numpy()[0]

n_CDS = candidates.select(
    (
        (pl.col('relation_to_gene') != 'intergenic') &
        pl.col('relation_to_gene').str.contains('protein_coding') &
        pl.col('relation_to_gene').str.contains('CDS')
    ).sum().alias('out')
)['out'].to_numpy()[0]

# non 3' and 5' UTR STRs
assert 0 == candidates.select(
    (
        (pl.col('relation_to_gene') != 'intergenic') &
        pl.col('relation_to_gene').str.contains('protein_coding') &
        ~pl.col('relation_to_gene').str.contains('CDS') &
        pl.col('relation_to_gene').str.contains('three_prime_UTR') &
        pl.col('relation_to_gene').str.contains('five_prime_UTR')
    ).sum().alias('out')
)['out'].to_numpy()[0]

# no unidentified UTR strs
assert 0 == candidates.select(
    (
        (pl.col('relation_to_gene') != 'intergenic') &
        pl.col('relation_to_gene').str.contains('protein_coding') &
        ~pl.col('relation_to_gene').str.contains('CDS') &
        pl.col('relation_to_gene').str.contains('UTR') &
        ~pl.col('relation_to_gene').str.contains('three_prime_UTR') &
        ~pl.col('relation_to_gene').str.contains('five_prime_UTR')
    ).sum().alias('out')
)['out'].to_numpy()[0]

n_five_UTR = candidates.select(
    (
        (pl.col('relation_to_gene') != 'intergenic') &
        pl.col('relation_to_gene').str.contains('protein_coding') &
        ~pl.col('relation_to_gene').str.contains('CDS') &
        pl.col('relation_to_gene').str.contains('five_prime_UTR')
    ).sum().alias('out')
)['out'].to_numpy()[0]

n_three_UTR = candidates.select(
    (
        (pl.col('relation_to_gene') != 'intergenic') &
        pl.col('relation_to_gene').str.contains('protein_coding') &
        ~pl.col('relation_to_gene').str.contains('CDS') &
        pl.col('relation_to_gene').str.contains('three_prime_UTR')
    ).sum().alias('out')
)['out'].to_numpy()[0]

# should be 1
n_exon_other = candidates.select(
    (
        (pl.col('relation_to_gene') != 'intergenic') &
        pl.col('relation_to_gene').str.contains('protein_coding') &
        ~pl.col('relation_to_gene').str.contains('CDS') &
        ~pl.col('relation_to_gene').str.contains('three_prime_UTR') &
        ~pl.col('relation_to_gene').str.contains('five_prime_UTR') &
        pl.col('relation_to_gene').str.contains('exon_other')
    ).sum().alias('out')
)['out'].to_numpy()[0]

# should be 1 more
n_intron = candidates.select(
    (
        (pl.col('relation_to_gene') != 'intergenic') &
        pl.col('relation_to_gene').str.contains('protein_coding') &
        ~pl.col('relation_to_gene').str.contains('CDS') &
        ~pl.col('relation_to_gene').str.contains('three_prime_UTR') &
        ~pl.col('relation_to_gene').str.contains('five_prime_UTR') &
        ~pl.col('relation_to_gene').str.contains('exon_other') &
        pl.col('relation_to_gene').str.contains('intron')
    ).sum().alias('out')
)['out'].to_numpy()[0]

'''
print(candidates.filter(
    (pl.col('relation_to_gene') != 'intergenic') &
    pl.col('relation_to_gene').str.contains('protein_coding') &
    ~pl.col('relation_to_gene').str.contains('CDS') &
    ~pl.col('relation_to_gene').str.contains('three_prime_UTR') &
    ~pl.col('relation_to_gene').str.contains('five_prime_UTR')
)['relation_to_gene'].to_numpy())

print('CDS', n_CDS, "3' UTR", n_three_UTR, "5' UTR", n_five_UTR, 'exon other', n_exon_other, 'n_coding_intronn', n_intron, 'noncoding', n_noncoding, 'intergenic', n_intergenic)
'''

# this isn't reliable if the data changes
categories=['Coding', "5' UTRs", "3' UTRs", "Other unannotated exon", "Intron", "lncRNA intron", "intergenic"]
histo = bokeh.plotting.figure(
    width=1200,
    height=400,
    title='Relationship between causal STR candidates and genes',
    x_range=categories,
    y_axis_label='count'
)

tops = [n_CDS, n_five_UTR, n_three_UTR, n_exon_other-1, n_intron+1, n_noncoding, n_intergenic]
cds = bokeh.models.ColumnDataSource(dict(
    cats = categories,
    idx=list(range(len(categories))),
    y=tops,
    text=[str(top) for top in tops]
))

histo.vbar(x='cats', top='y', width=0.9, source=cds)

'''
labels = bokeh.models.LabelSet(
    x=list(range(len(categories))), y=tops, text=[str(top) for top in tops],
    x_offset=5, y_offset=5, render_mode='canvas'
)
'''
labels = bokeh.models.LabelSet(
    x='idx', y='y', text='text', source=cds,
    x_offset=65, y_offset=5, level='glyph'
)
histo.add_layout(labels)
graphing_utils.export(histo, args.out, args.out.split('.')[-1])
