#!/usr/bin/env python3

import argparse
import ast
import os

import cyvcf2
import numpy as np
import pandas as pd
from statsmodels.regression.linear_model import WLS
from statsmodels.stats.weightstats import DescrStatsW

import python_array_utils as utils
import annotation_utils

ukb = os.environ['UKB']

def prep_dict(d):
    return {float(k): v for k, v in d.items() if not np.isnan(v)}

def rename_column_pd(df, old_name, new_name):
    df.rename(columns = {old_name: new_name}, inplace = True)

def main(outfname, readme, phenotype, literature_STRs, literature_STR_URL_list, cool_loci):
    with open(f'{ukb}/traits/phenotypes/{phenotype}_unit.txt') as unitfile:
        unit = next(unitfile).strip()

    col_dphen_unit = 'Δphenotype_per_additional_repeat_unit'
    col_dphen_sd = 'Δphenotype_per_s.d._increase_in_repeat_size'

    # ordered list of pairs (colname, description)
    columns = [
        ('chrom', None),
        ('signal_region', '{start bp}_{end bp} of region of association the STR resides in, '
         '1-based, inclusive'),
        ('start_pos', '1-based, inclusive'),
        ('end_pos', '1-based, inclusive'),
        ('SNPSTR_start_pos', 'the start position of the STR in the SNP-STR reference panel '
         '(may be smaller that start_pos if HipSTR called this STR with a physically phased '
         'SNP upstream of the STR)'),
        ('SNPSTR_ID', 'ID of the STR in the SNPSTR reference panel'),
        ('period', 'as specified in the SNPSTR reference panel'),
        ("repeat_unit", "as inferred by TRTools from the STR sequence and "
         "the period"),
        ('alleles', 'possible # copies of the repeat unit'),
        ('reference_allele', None),
        ('total_per_allele_dosages', 'sum of imputed dosage of each allele across '
         'both chromsomes in all samples'),
        ('total_hardcall_alleles', 'number of each allele in the population using '
         "Beagle's called most probable phased genotypes"),
        ('total_hardcall_genotypes', 'number of each unphased genotype in the '
         "population using Beagle's called most probable phased genotypes"),
        ('subset_total_per_allele_dosages', 'as total_per_allele_dosages, but '
         'calculated on the subset of samples used in the association'),
        ('subset_total_hardcall_alleles', 'as total_hardcall_alleles , but '
         'calculated on the subset of samples used in the association'),
        ('subset_total_hardcall_genotypes', 'as total_hardcall_genotypes, but '
         'calculated on the subset of samples used in the association'),
        ('subset_multiallelicness', '% of dosage attributed to alleles that are not '
         'the first or second most frequent alleles, by dosage, within the subset of '
         'samples used in the association'),
        ('subset_heterozygosity', 'the expected percent of heterozygous samples in the '
         'association sample subset as calculated from allele frequencies'),
        ('subset_entropy', 'the entropy in the association sample subset calculated '
         'from allele frequencies'),
        ('subset_HWEP', 'the Hardy-Weinberg p-value in the association sample subset '
         'comparing hardcall genotypes and allele frequencies'),
        ('subset_allele_dosage_r2', 'a metric of imputation accuracy, the per-allele R^2 '
         'between dosage-weighted allele length and allele hardcalls across samples'),
        ('association_p_value', 'linear regression of rank inverse normalized phenotype '
         "values vs repeat length on QC'ed sample subset, with covariates"),
        (col_dphen_unit, 'linear regression on raw phenotypes vs repeat length on '
         f"QC'ed sample subset, no covariates included, phenotype measured in {unit}"),
        (col_dphen_sd, 'linear regression on raw phenotypes vs repeat length on '
         f"QC'ed sample subset, no covariates included, phenotype measured in {unit}"),
        ('pcausal', 'FINEMAP posterior probability of causality'),
        ('mentioned_in_literature', 'Whether or not this we know of a citation '
         'saying this STR is likely causal for this trait. False here means we have '
         'not checked whether or not this is reported in the literature, not that it has not.'),
        ('literature_inclusion_url', 'If mentioned_in_literature, then the corresponding ULR. '
          'Otherwise NA. '),
        ('included_only_due_to_literature', 'NA if not mentioned in literature. True if the locus '
         'would have been filtered, but is included due to literature. False otherwise. '),
        ('curated', 'Whether or not manual examination was used to decide this was an '
         'exciting locus worthy of extra attention. '),
        ('included_only_due_to_curation', 'NA if not curated. True if the locus '
         'would have been filtered, but is included due to manual curation. '
         'False if the locus was curated but would have been included regardless (not filtered). '),
        ('nearby_exons', "A comma separated list of distance:'upstream'/'downstream'/'inside':"
         "exon-start:exon-end:gene-name:gene-type, where upstream means upstream of "
         "the exon in that gene's direction. All exons within a 1000bp radius, or 'none'. (See "
         'https://www.gencodegenes.org/pages/biotypes.html for gene type meanings.) '),
        ('closest_gene', "distance:'upstream'/'downstream'/'inside':gene-name:gene-type. Possibly a comma "
         'separated list if multiple genes are tied by distance. '
         "(See https://www.gencodegenes.org/pages/biotypes.html for gene type meanings.)"),
        ('nearby_genes', 'A comma separated list of entries like that in closest_gene. '
         "All genes within a 100kbp radius, or 'none'."),
        ('relation_to_gene', "If this STR intersects a single gene with a single category, then "
         "'CDS'/'5_prime_UTR'/'3_prime_UTR'/'unannotated_UTR'/'exon_other'/'intron':gene-name:gene-type. "
         "If it intersects multiple categories of the same gene, then those categories are "
         "comma separated before the gene and gene type. If it intersects multiple genes, then "
         "'multigene;' followed by the categorization for each of the genes as above, separated "
         "by semicolons. If it intersects no genes, then 'intergenic'. If it partially intersects "
         "at least one gene, then 'intergenic;' followed by the information above. "
         "Note that STRs which have "
         "only exonic annotations might still be introns in relation to other transcripts of the "
         "gene they overlap. "
         "(See https://www.gencodegenes.org/pages/biotypes.html for gene type meanings.) "),
        ('transcribed', "Does this intersect a transcript? If so, 'transcribed;' followed by "
         "a comma separated list of transcript-name:transcript-type:transcript-support-level. "
         "Otherwise 'untranscribed'."
         "(See https://www.gencodegenes.org/pages/biotypes.html for transcript-type meanings "
         "and https://www.gencodegenes.org/pages/data_format.html for transcript-support-level "
         "meanings.)"),
    ]

    readme.write(
        'table.tab contains one row for each variant that was not filtered prior to association '
        'testing and either was reported as an association in literature, or I manually marked it as interesting '
        'or both had association p-value >= 5e-8 and '
        'FINEMAP posterior probability of causaility >= 0.05\n'
        '\n'
        'table.tab contains the following columns:\n'
    )
    for column, description in columns:
        readme.write(column)
        if description is not None:
            readme.write(' - ')
            readme.write(description)
        readme.write('\n')
    readme.write('\n')

    print("Reading finemapping output ...", flush=True)
    # chrom, signal_region, SNPSTR_start_pos, pcausal
    signals_fname = f'{ukb}/finemapping/finemap_results/{phenotype}/summary/all_STR_contribs.tab'
    signals = pd.read_csv(
        signals_fname,
        skiprows = 1,
        delimiter='\t',
        dtype=utils.get_dtypes(signals_fname)
    )
    rename_column_pd(signals, 'signal', 'signal_region')

    signal_split = signals['signal_region'].str.split('_', n=1, expand=True)
    signals['chrom'] = signal_split[0].astype(int)
    str_split = signals['STR'].str.split('_', n=1, expand=True)
    signals['pos'] = str_split[1].astype(int)

    str_results_fname = f'{ukb}/association/results/{phenotype}/my_str/results.tab'

    '''
    chrom
    SNPSTR_start_pos
    repeat_unit
    period
    alleles
    ref_len
    total_per_allele_dosages
    total_hardcall_alleles
    total_hardcall_genotypes
    subset_total_per_allele_dosages
    subset_total_hardcall_alleles
    subset_total_hardcall_genotypes
    subset_het
    subset_entropy
    subset_HWEP
    subset_allele_dosage_r2
    association_p_value
    '''
    print("Reading associations ...", flush=True)
    associations = pd.read_csv(
        str_results_fname,
        header=0,
        delimiter='\t',
        dtype=utils.get_dtypes(str_results_fname, {'locus_filtered': str})
    )
    rename_column_pd(associations, 'ref_len', 'reference_allele')
    rename_column_pd(associations, 'motif', 'repeat_unit')
    rename_column_pd(associations, 'subset_het', 'subset_heterozygosity')
    rename_column_pd(associations, f'p_{phenotype}', 'association_p_value')
    for STR in literature_STRs:
        chrom, pos = STR.split(':')
        if not np.any((associations['chrom'] == int(chrom)) & (associations['pos'] == int(pos))):
            readme.write(f'Never ran association on literature reported STR association at {STR}\n')
            literature_STRs[STR] = True

    print("Merging ...", flush=True)
    signals = signals.merge(
        associations,
        on=['chrom', 'pos'],
        how='right'
    )
    literature_STR_rows = None
    literature_STR_URLs = ['NA']*len(signals)
    for (STR, already_dropped), URL in zip(literature_STRs.items(), literature_STR_URL_list):
        if already_dropped:
            continue
        chrom, pos = STR.split(':')
        new_rows = (signals['chrom'] == int(chrom)) & (signals['pos'] == int(pos))
        assert np.sum(new_rows) == 1
        for idx in np.nonzero(new_rows.to_numpy())[0]:
            literature_STR_URLs[idx] = URL
        if literature_STR_rows is None:
            literature_STR_rows = new_rows
        else:
            literature_STR_rows |= new_rows
    if literature_STR_rows is None:
        literature_STR_rows = np.full((signals.shape[0]), False, dtype=bool)
    else:
        literature_STR_rows = literature_STR_rows.to_numpy()

    signals['mentioned_in_literature'] = literature_STR_rows
    signals['literature_inclusion_url'] = literature_STR_URLs
    only_lit = ['NA']*signals.shape[0]
    for idx in np.nonzero(literature_STR_rows)[0]:
        only_lit[idx] = 'True'
    for idx in np.nonzero(
        literature_STR_rows &
        ~np.isnan(signals['pcausal'].to_numpy()) &
        (signals['pcausal'].to_numpy() >= 0.05)
    )[0]:
        only_lit[idx] = 'False'

    signals['included_only_due_to_literature'] = only_lit

    cool_loci_rows = None
    for STR in cool_loci:
        chrom, pos = STR.split(':')
        new_rows = (signals['chrom'] == int(chrom)) & (signals['pos'] == int(pos))
        assert np.sum(new_rows) == 1
        if cool_loci_rows is None:
            cool_loci_rows = new_rows
        else:
            cool_loci_rows |= new_rows
    if cool_loci_rows is None:
        cool_loci_rows = np.full((signals.shape[0]), False, dtype=bool)
    else:
        cool_loci_rows = cool_loci_rows.to_numpy()

    signals['curated'] = cool_loci_rows
    only_cool = ['NA']*signals.shape[0]
    for idx in np.nonzero(cool_loci_rows)[0]:
        only_cool[idx] = 'True'
    for idx in np.nonzero(
        cool_loci_rows &
        ~np.isnan(signals['pcausal'].to_numpy()) &
        (signals['pcausal'].to_numpy() >= 0.05)
    )[0]:
        only_cool[idx] = 'False'

    signals['included_only_due_to_curation'] = only_cool

    signals = signals[
        literature_STR_rows |
        cool_loci_rows | (
            ~np.isnan(signals['pcausal'].to_numpy()) &
            (signals['pcausal'].to_numpy() >= 0.05) &
            (signals['association_p_value'].to_numpy() < 5e-8)
    )]

    signals.reset_index(inplace=True)
    nrows = signals.shape[0]

    rename_column_pd(signals, 'pos', 'SNPSTR_start_pos')

    # calculate effect sizes and third on allele freqs
    print("Working on effect sizes and multiallelicness ...", flush=True)
    dphen_units = np.full(nrows, np.nan)
    dphen_sds = np.full(nrows, np.nan)
    third_on_allele_freqs = np.full(nrows, np.nan)
    for idx in range(nrows):
        single_dosages = prep_dict(ast.literal_eval(signals['subset_total_per_allele_dosages'][idx]))
        single_dosages = prep_dict(ast.literal_eval(signals['subset_total_per_allele_dosages'][idx]))
        total_dosage = sum(single_dosages.values())
        single_freqs = {key: value/total_dosage for key, value in single_dosages.items()}

        if len(single_freqs) <= 2:
            third_on_allele_freqs[idx] = 0
        else:
            third_on_allele_freqs[idx] = sum(sorted(single_freqs.values())[:-2])

        doubled_freqs = {}
        for k1 in single_freqs:
            for k2 in single_freqs:
                key = k1 + k2
                if key not in doubled_freqs:
                    doubled_freqs[key] = 0
                doubled_freqs[key] += single_freqs[k1]*single_freqs[k2]

        phenotype_vals = prep_dict(ast.literal_eval(signals[f'mean_{phenotype}_per_single_dosage'][idx]))
        for key in phenotype_vals:
            if key not in doubled_freqs:
                print(phenotype_vals)
                print(doubled_freqs)
                assert False
        for key in doubled_freqs:
            if key not in phenotype_vals:
                if not doubled_freqs[key] <= 0.00001:
                    print(phenotype_vals)
                    print(doubled_freqs)
                    assert False

        xs = sorted(phenotype_vals.keys())
        ys = np.array([phenotype_vals[x] for x in xs])
        weights = np.array([doubled_freqs[x] for x in xs])
        xs = np.array(xs).reshape(-1, 1)
        xs = np.concatenate((xs, np.ones((xs.shape[0], 1))), axis=1)

        model = WLS(ys, xs, weights)
        reg_result = model.fit()
        coef = reg_result.params[0]
        scaled_coef = coef*DescrStatsW(xs[:, 0], weights=weights).std
        dphen_units[idx] = coef
        dphen_sds[idx] = scaled_coef
    signals['subset_multiallelicness'] = third_on_allele_freqs
    signals[col_dphen_unit] = dphen_units
    signals[col_dphen_sd] = dphen_sds

    print("Reading info vcfs ...", flush=True)
    chroms = []
    SNPSTR_poses = []
    SNPSTR_IDs = []
    starts = []
    ends = []
    for chrom in range(1, 23):
        vcf = cyvcf2.VCF(f'{ukb}/snpstr/info_field/chr{chrom}.vcf.gz')
        for variant in vcf:
            chroms.append(int(variant.CHROM))
            SNPSTR_poses.append(variant.POS)
            SNPSTR_IDs.append(variant.ID)
            starts.append(variant.INFO['START'])
            ends.append(variant.INFO['END'])
        vcf.close()

    info_df = pd.DataFrame.from_dict({
        'chrom': chroms,
        'SNPSTR_start_pos': SNPSTR_poses,
        'SNPSTR_ID': SNPSTR_IDs,
        'start_pos': starts,
        'end_pos': ends
    })

    signals = signals.merge(
        info_df,
        on=['chrom', 'SNPSTR_start_pos'],
        how="left"
    )

    # handle gene & transcript annotations
    print("Loading gene and transcript annotations ...", flush=True)
    signals['pos'] = signals['start_pos']
    nearby_exon_merge = annotation_utils.get_merged_annotations(signals, f'{ukb}/side_analyses/str_annotations/nearby_exon_d_1000')
    closest_gene_merge = annotation_utils.get_merged_annotations(signals, f'{ukb}/side_analyses/str_annotations/closest_gene', distance=True)
    nearby_gene_merge = annotation_utils.get_merged_annotations(signals, f'{ukb}/side_analyses/str_annotations/nearby_gene_d_100000')

    print("Loading high level intersections...", flush=True)
    transcript_intersect_merge = annotation_utils.get_merged_annotations(signals, f'{ukb}/side_analyses/str_annotations/intersects_transcript', bp_overlap=True)
    gene_intersect_merge = annotation_utils.get_merged_annotations(signals, f'{ukb}/side_analyses/str_annotations/intersects_gene', bp_overlap=True)
    exon_intersect_merge = annotation_utils.get_merged_annotations(signals, f'{ukb}/side_analyses/str_annotations/intersects_exon', bp_overlap=True)

    print("Loading low level intersections...", flush=True)
    CDS_intersect_merge = annotation_utils.get_merged_annotations(signals, f'{ukb}/side_analyses/str_annotations/intersects_CDS', bp_overlap=True)
    five_prime_UTR_intersect_merge = annotation_utils.get_merged_annotations(signals, f'{ukb}/side_analyses/str_annotations/intersects_five_prime_UTR', bp_overlap=True)
    three_prime_UTR_intersect_merge = annotation_utils.get_merged_annotations(signals, f'{ukb}/side_analyses/str_annotations/intersects_three_prime_UTR', bp_overlap=True)
    UTR_intersect_merge = annotation_utils.get_merged_annotations(signals, f'{ukb}/side_analyses/str_annotations/intersects_UTR', bp_overlap=True)

    sub_exon_types = pd.concat([CDS_intersect_merge, five_prime_UTR_intersect_merge, three_prime_UTR_intersect_merge, UTR_intersect_merge])

    nrows = signals.shape[0]

    nearby_exons = ['none']*nrows
    closest_gene = [None]*nrows
    nearby_genes = ['none']*nrows
    relation_to_gene = ['']*nrows
    transcribed = ['untranscribed']*nrows

    print("Handling gene and transcript annotations ...", flush=True)
    for idx in range(nrows):
        chrom = signals['chrom'][idx]
        snpstr_pos = signals['SNPSTR_start_pos'][idx]
        start_pos = signals['start_pos'][idx]
        end_pos = signals['end_pos'][idx]
        for line in nearby_exon_merge[
            (nearby_exon_merge['chrom'] == chrom) & (nearby_exon_merge['STR_pos'] == snpstr_pos)
        ].itertuples():
            if nearby_exons[idx] == 'none':
                nearby_exons[idx] = ''
            else:
                nearby_exons[idx] += ','
            nearby_exons[idx] += annotation_utils.distance(
                start_pos, end_pos, line.annotation_pos, line.annotation_end_pos
            ) + ":"
            nearby_exons[idx] += annotation_utils.get_relation(
                start_pos, end_pos, line.annotation_pos, line.annotation_end_pos, line.annotation_strand
            ) + ":"
            nearby_exons[idx] += str(line.annotation_pos) + ":"
            nearby_exons[idx] += str(line.annotation_end_pos) + ":"
            nearby_exons[idx] += annotation_utils.get_gff_kvp(line.annotation_info, 'gene_name') + ":"
            nearby_exons[idx] += annotation_utils.get_gff_kvp(line.annotation_info, 'gene_type')

        for line in closest_gene_merge[
            (closest_gene_merge['chrom'] == chrom) & (closest_gene_merge['STR_pos'] == snpstr_pos)
        ].itertuples():
            if closest_gene[idx] is not None:
                closest_gene[idx] += ','
            else:
                closest_gene[idx] = ''
            closest_gene[idx] += str(line.annotation_distance) + ":"
            closest_gene[idx] += annotation_utils.get_relation(
                start_pos, end_pos, line.annotation_pos, line.annotation_end_pos, line.annotation_strand
            ) + ":"

            closest_gene[idx] += annotation_utils.get_gff_kvp(line.annotation_info, 'gene_name') + ":"
            closest_gene[idx] += annotation_utils.get_gff_kvp(line.annotation_info, 'gene_type')

        for line in nearby_gene_merge[
            (nearby_gene_merge['chrom'] == chrom) & (nearby_gene_merge['STR_pos'] == snpstr_pos)
        ].itertuples():
            if nearby_genes[idx] == 'none':
                nearby_genes[idx] = ''
            else:
                nearby_genes[idx] += ','
            nearby_genes[idx] += annotation_utils.distance(
                start_pos, end_pos, line.annotation_pos, line.annotation_end_pos
            ) + ":"
            nearby_genes[idx] += annotation_utils.get_relation(
                start_pos, end_pos, line.annotation_pos, line.annotation_end_pos, line.annotation_strand
            ) + ":"
            nearby_genes[idx] += annotation_utils.get_gff_kvp(line.annotation_info, 'gene_name') + ":"
            nearby_genes[idx] += annotation_utils.get_gff_kvp(line.annotation_info, 'gene_type')

        # dict from gene to types of intersections
        gene_intersections = {}
        partial_overlap = False
        for line in gene_intersect_merge[
            (gene_intersect_merge['chrom'] == chrom) & (gene_intersect_merge['STR_pos'] == snpstr_pos)
        ].itertuples():
            gene_name = annotation_utils.get_gff_kvp(line.annotation_info, 'gene_name')
            gene_type = annotation_utils.get_gff_kvp(line.annotation_info, 'gene_type')
            if gene_name in gene_intersections:
                # already described this gene, just hope the description is the same
                continue
            gene_intersections[gene_name] = [
                gene_type,
                {p: False for p in range(start_pos, end_pos+1)}, # bps covered by gene annotation
                {p: False for p in range(start_pos, end_pos+1)}, # bps covered by exon annotations
                set(), # types of sub-exon intersections
                {p: False for p in range(start_pos, end_pos+1)} # bps covered by sub-exon annotations
            ]
            for p in gene_intersections[gene_name][1]:
                if line.annotation_pos <= p <= line.annotation_end_pos:
                    gene_intersections[gene_name][1][p] = True
                else:
                    partial_overlap = True

        if len(gene_intersections) == 0:
            relation_to_gene[idx] = 'intergenic'

        if partial_overlap:
            relation_to_gene[idx] = 'intergenic;'

        if len(gene_intersections) > 1:
            relation_to_gene[idx] += 'multigene;'

        for line in exon_intersect_merge[
            (exon_intersect_merge['chrom'] == chrom) & (exon_intersect_merge['STR_pos'] == snpstr_pos)
        ].itertuples():
            gene_name = annotation_utils.get_gff_kvp(line.annotation_info, 'gene_name')
            if gene_name not in gene_intersections:
                raise ValueError(f"STR {chrom} {start_pos} gene name {gene_name} not in list of overlapping genes {set(gene_intersections.keys())}")
            for p in gene_intersections[gene_name][2]:
                if line.annotation_pos <= p <= line.annotation_end_pos:
                    gene_intersections[gene_name][2][p] = True

        for line in sub_exon_types[
            (sub_exon_types['chrom'] == chrom) & (sub_exon_types['STR_pos'] == snpstr_pos)
        ].itertuples():
            gene_name = annotation_utils.get_gff_kvp(line.annotation_info, 'gene_name')
            if gene_name not in gene_intersections:
                raise ValueError(f"STR {chrom} {start_pos} gene name {gene_name} not in list of overlapping genes {set(gene_intersections.keys())}")
            if gene_intersections[gene_name][2] == False:
                raise ValueError(f"STR {chrom} {start_pos} intersects an exon type but isn't in an exon!")
            gene_intersections[gene_name][3].add(line.annotation_feature_type)
            for p in gene_intersections[gene_name][4]:
                if line.annotation_pos <= p <= line.annotation_end_pos:
                    if not gene_intersections[gene_name][2][p]:
                        raise ValueError(f"STR {chrom} {start_pos} at position {p} has sub-exon type {line.annotation_feature_type} but is not exonic!")
                    gene_intersections[gene_name][4][p] = True

        first = True
        for gene_name, gene_data in gene_intersections.items():
            if not first:
                relation_to_gene[idx] += ';'
            first = False
            types = gene_data[3]
            for p, is_gene in gene_data[1].items():
                if is_gene and not gene_data[2][p]:
                    types.add('intron')
                    break
            for p, is_gene in gene_data[2].items():
                if is_gene and not gene_data[4][p]:
                    types.add('exon_other')
                    break
            relation_to_gene[idx] += ','.join(types) + ":"
            relation_to_gene[idx] += gene_name + ":" + gene_data[0]

        for line in transcript_intersect_merge[
            (transcript_intersect_merge['chrom'] == chrom) & (transcript_intersect_merge['STR_pos'] == snpstr_pos)
        ].itertuples():
            if transcribed[idx] == 'untranscribed':
                transcribed[idx] = 'transcribed;'
            else:
                transcribed[idx] += ','
            transcribed[idx] += annotation_utils.get_gff_kvp(line.annotation_info, 'transcript_name') + ":"
            transcribed[idx] += annotation_utils.get_gff_kvp(line.annotation_info, 'transcript_type') + ":"
            transcribed[idx] += annotation_utils.get_gff_kvp(line.annotation_info, 'transcript_support_level')

    signals['nearby_exons'] = nearby_exons
    signals['closest_gene'] = closest_gene
    signals['nearby_genes'] = nearby_genes
    signals['relation_to_gene'] = relation_to_gene
    signals['transcribed'] = transcribed

    signals = signals[[colname for colname, description in columns]]

    signals.to_csv(
        f'{outfname}.tab',
        sep='\t',
        index=False,
        na_rep = 'NA'
    )

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('phenotype')
    parser.add_argument('outfname') # writes out ${outfname}.tab and ${outfname}_README.tab
    parser.add_argument('--previous-STR-findings', nargs='+', default=[])
    parser.add_argument('--previous-STR-finding-URLs', nargs='+', default=[])
    parser.add_argument('--cool-loci', nargs='+', default=[])
    args = parser.parse_args()
    phenotype = args.phenotype
    outfname = args.outfname
    previous_STRs = {STR : False for STR in args.previous_STR_findings}
    previous_STR_URLs = args.previous_STR_finding_URLs

    assert len(previous_STRs) == len(previous_STR_URLs)

    with open(f'{outfname}_README.txt', 'w') as readme:
        main(outfname, readme, phenotype, previous_STRs, previous_STR_URLs, args.cool_loci)
