#!/usr/bin/env python3

import argparse
import datetime

import polars as pl

import sample_utils

def write_input_variants(workdir, outdir, gts_dir, plink_results_fname, str_results_fname, filter_set_fname, samples_fname, readme, phenotype, chrom, start, end, inclusion_threshold, mac, snp_str_ratio, total_prob):
    '''
    write README.txt
    write finemap_input.z
    write finemap_input.master
    '''

    n_samples = sample_utils.n_samples(samples_fname)

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
    if str_results_fname != "None":
        strs = pl.scan_csv(
            str_results_fname, sep='\t', dtypes={'locus_filtered': str}
        ).filter(
            (pl.col('chrom') == chrom) &
            (pl.col('pos') >= start) &
            (pl.col('pos') <= end) &
            (pl.col('locus_filtered') == 'False') &
            (pl.col(f'p_{phenotype}') < inclusion_threshold)
        ).filter(
            # STRs with duplicate loci that shouldn't have been in the reference panel
            ((pl.col('chrom') != 17) | (pl.col('pos') != 80520458)) &
            ((pl.col('chrom') != 1) | (pl.col('pos') != 247747217)) &
            ((pl.col('chrom') != 1) | (pl.col('pos') != 247848392)) &
            ((pl.col('chrom') != 21) | (pl.col('pos') != 47741815)) &
            ((pl.col('chrom') != 8) | (pl.col('pos') != 145231731))
        ).collect()
        if strs.shape[0] > 0:
            strs = strs.select([
                ('STR_' + pl.col('pos').cast(str)).alias('rsid'),
                ('0' + pl.col('chrom').cast(str)).str.slice(-2).alias('chromosome'),
                pl.col('pos').alias('position'),
                pl.lit('nan').alias('allele1'),
                pl.lit('nan').alias('allele2'),
                pl.lit('nan').alias('maf'),
                pl.col(f'coeff_{phenotype}').alias('beta'),
                pl.col(f'se_{phenotype}').alias('se')
            ])
        else:
            strs = strs.select([
                pl.lit('empty').alias('rsid'),
                pl.lit('empty').alias('chromosome'),
                pl.col('pos').alias('position'),
                pl.lit('nan').alias('allele1'),
                pl.lit('nan').alias('allele2'),
                pl.lit('nan').alias('maf'),
                pl.col(f'coeff_{phenotype}').alias('beta'),
                pl.col(f'se_{phenotype}').alias('se')
            ])

        if mac:
            strs = strs.filter(~pl.col('position').is_in(strs_exclude_mac))

        assert strs.unique(subset=['chromosome', 'position']).shape[0] == strs.shape[0]

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

    if str_results_fname != "None":
        vars_df = pl.concat([strs, snps])
    else:
        vars_df = snps

    if total_prob is not None:
        vars_df = vars_df.with_column(
            pl.lit(total_prob/(n_snps + n_strs)).alias('prob')
        )

    vars_df.write_csv(f'{workdir}/finemap_input.z', sep=' ')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('outdir')
    parser.add_argument('gts_dir')
    parser.add_argument('plink_results')
    parser.add_argument('str_results')
    parser.add_argument('variants_to_filter')
    parser.add_argument('samples_fname')
    parser.add_argument('phenotype')
    parser.add_argument('chrom', type=int)
    parser.add_argument('start_pos', type=int)
    parser.add_argument('end_pos', type=int)
    parser.add_argument('--snp-str-ratio', type=float, default=None)
    parser.add_argument('--total-prob', type=float, default=None)
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
    with open(f'{outdir}/README.txt', 'w') as readme:
        write_input_variants(outdir, outdir, gts_dir, args.plink_results, args.str_results, args.variants_to_filter, args.samples_fname, readme, phenotype, chrom, start_pos, end_pos, args.inclusion_threshold, args.mac, args.snp_str_ratio, args.total_prob)

if __name__ == '__main__':
    main()

