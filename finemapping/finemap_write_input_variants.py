#!/usr/bin/env python3

import argparse
import datetime
import os
import tempfile

import numpy as np
import polars as pl

import python_file_utils as file_utils
import sample_utils

ukb = os.environ['UKB']

def write_input_variants(workdir, outdir, gts_dir, readme, phenotype, chrom, start, end, inclusion_threshold, mac, snp_str_ratio, total_prob, use_PACSIN2):
    '''
    write README.txt
    write finemap_input.z
    write finemap_innput.master
    '''

    sample_idx = sample_utils.get_samples_idx_phenotype('white_brits', phenotype)
    n_samples = np.sum(sample_idx)

    if mac:
        mac_threshold = int(mac[0])
        snp_mac_fname = mac[1]
        str_mac_fname = mac[2]
        snps_exclude_mac = pl.scan_csv(
            snp_mac_fname,
            sep='\t'
        ).filter(
            pl.col('ALT_CTS') < mac_threshold
        ).select(
            ('SNP_' + pl.col('#POS').cast(str) + '_' + pl.col('REF') + '_' + pl.col('ALT')).alias('varname')
        ).collect()['varname'].to_list()
        # need to make that look like a list of strings to polars b/c buggy, so add a single nonsense to it
        snps_exclude_mac.append('asdf')

        strs_exclude_mac = pl.scan_csv(
            str_mac_fname,
            sep='\t'
        ).filter(
            pl.col('mac') < mac_threshold
        ).select(
            'pos'
        ).collect()['pos'].to_list()

    plink_results_fname = f'{ukb}/association/results/{phenotype}/plink_snp/results.tab'
    str_results_fname = f'{ukb}/association/results/{phenotype}/my_str/results.tab'
    filter_set_fname = f'{ukb}/finemapping/str_imp_snp_overlaps/chr{chrom}_to_filter.tab'

    with open(f'{workdir}/finemap_input.master', 'w') as finemap_master:
        finemap_master.write(
            'z;ld;snp;config;cred;log;n_samples\n'
            f'{outdir}/finemap_input.z;'
            f'{gts_dir}/all_variants.ld;'
            f'{outdir}/finemap_output.snp;'
            f'{outdir}/finemap_output.config;'
            f'{outdir}/finemap_output.cred;'
            f'{outdir}/finemap_output.log;'
            f'{n_samples}'
        )

    today = datetime.datetime.now().strftime("%Y_%M_%D")
    readme.write(
        f'Run date: {today}\n'
        'Manually generating variant-variant LD for each imputed SNP each STR in the region '
        'where an association was successfully '
        f'performed and had p < {inclusion_threshold} and the SNP was not in the filter set\n'
        f'(Filter set at {filter_set_fname})\n'
        'Correlation is STR length dosage vs SNP dosage.\n'
        'Running FINEMAP with that list of imputed SNPs and STRs.\n'
    )

    # load STRs
    strs = pl.scan_csv(
        str_results_fname, sep='\t', dtypes={'locus_filtered': str}
    ).filter(
        (pl.col('chrom') == chrom) &
        (pl.col('pos') >= start) &
        (pl.col('pos') <= end) &
        (pl.col('locus_filtered') == 'False') &
        (pl.col(f'p_{phenotype}') < inclusion_threshold)
    ).select([
        ('STR_' + pl.col('pos').cast(str)).alias('rsid'),
        ('0' + pl.col('chrom').cast(str)).str.slice(-2).alias('chromosome'),
        pl.col('pos').alias('position'),
        pl.lit('nan').alias('allele1'),
        pl.lit('nan').alias('allele2'),
        pl.lit('nan').alias('maf'),
        pl.col(f'coeff_{phenotype}').alias('beta'),
        pl.col(f'se_{phenotype}').alias('se'),
    ]).collect()

    if mac:
        strs = strs.filter(~pl.col('position').is_in(strs_exclude_mac))

    if use_PACSIN2:
        strs = strs.filter(
            pl.col('pos') != 43385872
        )
        pacsin2_strs = pl.read_csv(
            f'{ukb}/association/spot_test/white_brits/{phenotype}/PACSIN2.tab',
            sep='\t'
        ).filter(
            pl.col('pos').is_in([43385866, 43385875, 43385893])
        ).select([
            ('PACSIN2_STR_' + pl.col('pos').cast(str)).alias('rsid'),
            ('0' + pl.col('chrom').cast(str)).str.slice(-2).alias('chromosome'),
            pl.col('pos').alias('position'),
            pl.lit('nan').alias('allele1'),
            pl.lit('nan').alias('allele2'),
            pl.lit('nan').alias('maf'),
            pl.col(f'coeff_{phenotype}').alias('beta'),
            pl.col(f'se_{phenotype}').alias('se'),
        ])
        strs = pl.concat([strs, pacsin2_strs])

    assert strs.shape[0] > 0
    assert strs.distinct(subset=['chromosome', 'position']).shape[0] == strs.shape[0]

    n_strs = strs.shape[0]

    # load SNPs
    snps_to_filter = set()
    with open(filter_set_fname) as filter_file:
        next(filter_file) # skip header
        for line in filter_file:
            pos, ref, alt = line.strip().split('\t')[3:6]
            snps_to_filter.add(f'{pos}_{ref}_{alt}')

    snps = pl.scan_csv(
        plink_results_fname, sep='\t', null_values='NA'
    ).filter(
        (pl.col('#CHROM') == chrom) &
        (pl.col('POS') >= start) &
        (pl.col('POS') <= end) &
        (pl.col('ERRCODE') == '.') &
        (pl.col('P') < inclusion_threshold) &
        ~(pl.col('POS').cast(str) + '_' + pl.col('REF') + '_' + pl.col('ALT')).is_in(list(snps_to_filter))
    ).select([
        ('SNP_' + pl.col('POS').cast(str) + '_' + pl.col('REF') + '_' + pl.col('ALT')).alias('rsid'),
        ('0' + pl.col('#CHROM').cast(str)).str.slice(-2).alias('chromosome'),
        pl.col('POS').alias('position'),
        pl.col('REF').alias('allele1'),
        pl.col('ALT').alias('allele2'),
        pl.lit('nan').alias('maf'),
        pl.col('BETA').alias('beta'),
        pl.col('SE').alias('se'),
    ]).collect()

    if mac:
        snps=snps.filter(~pl.col('rsid').is_in(snps_exclude_mac))

    n_snps = snps.shape[0]

    if snp_str_ratio is not None:
        strs = strs.with_column(
            pl.lit(1/(n_strs + snp_str_ratio*n_snps)).alias('prob')
        )
        snps = snps.with_column(
            pl.lit(snp_str_ratio/(n_strs + snp_str_ratio*n_snps)).alias('prob')
        )

    vars_df = pl.concat([strs, snps])

    if total_prob is not None:
        vars_df = vars_df.with_column(
            pl.lit(total_prob/(n_snps + n_strs)).alias('prob')
        )

    vars_df.to_csv(f'{workdir}/finemap_input.z', sep=' ')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('outdir')
    parser.add_argument('gts_dir')
    parser.add_argument('phenotype')
    parser.add_argument('chrom', type=int)
    parser.add_argument('start_pos', type=int)
    parser.add_argument('end_pos', type=int)
    parser.add_argument('--snp-str-ratio', type=float, default=None)
    parser.add_argument('--total-prob', type=float, default=None)
    parser.add_argument('--three-PACSIN2-STRs', action='store_true', default=False)
    parser.add_argument('--inclusion-threshold', type=float, default=0.05)
    parser.add_argument('--mac', nargs=3, default=None)
    args = parser.parse_args()

    assert (args.snp_str_ratio is not None) + (args.total_prob is not None) <= 1

    phenotype = args.phenotype
    chrom = args.chrom
    start_pos = args.start_pos
    end_pos = args.end_pos
    assert start_pos < end_pos

    outdir = args.outdir
    gts_dir= args.gts_dir
    with file_utils.temp_dir(outdir) as tempdir:
        with open(f'{tempdir}/README.txt', 'w') as readme:
            write_input_variants(tempdir, outdir, gts_dir, readme, phenotype, chrom, start_pos, end_pos, args.inclusion_threshold, args.mac, args.snp_str_ratio, args.total_prob, args.three_PACSIN2_STRs)
        file_utils.move_files(tempdir, outdir)

if __name__ == '__main__':
    main()

