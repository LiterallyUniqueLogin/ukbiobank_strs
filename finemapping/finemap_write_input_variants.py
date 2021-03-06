#!/usr/bin/env python3

import argparse
import csv
import datetime
import os
import pathlib

import python_file_utils as file_utils

ukb = os.environ['UKB']

inclusion_threshold = 0.05

def write_input_variants(workdir, readme, phenotype, chrom, start_pos, end_pos):
    '''
    write README.txt
    write finemap_input.z
    sometimes write no_strs
    '''
    plink_results_fname = f'{ukb}/association/results/{phenotype}/plink_snp/results.tab'
    str_results_fname = f'{ukb}/association/results/{phenotype}/my_str/results.tab'
    filter_set_fname = f'{ukb}/finemapping/str_imp_snp_overlaps/chr{chrom}_to_filter.tab'

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

    with open(f'{workdir}/finemap_input.z', 'w') as finemap_input_z:
        finemap_input_z.write('rsid chromosome position allele1 allele2 maf beta se\n')

        any_strs = False
        prev_str_pos = None
        # load STRs
        with open(str_results_fname) as str_results_file:
            str_results_reader = csv.reader(str_results_file, delimiter='\t')
            header = next(str_results_reader)
            cols = {
                col: header.index(col) for col in
                ['chrom', 'pos', 'locus_filtered', f'p_{phenotype}', f'coeff_{phenotype}', f'se_{phenotype}']
            }

            # assumes ordered numeric chromosomes
            for result in str_results_reader:
                result_chrom = int(result[cols['chrom']])
                result_pos = int(result[cols['pos']])
                if (result_chrom, result_pos) < (chrom, start_pos):
                    continue
                if (result_chrom, result_pos) > (chrom, end_pos):
                    break
                if result[cols['locus_filtered']] != 'False' or float(result[cols[f'p_{phenotype}']]) >= inclusion_threshold:
                    continue
                if result_pos == prev_str_pos:
                    raise ValueError(f"Two STR poses at the same location {result_pos}!")
                prev_str_pos = result_pos
                any_strs = True

                beta = result[cols[f'coeff_{phenotype}']]
                se = result[cols[f'se_{phenotype}']]
                # I have forced there to be a unique STR per position during association testing
                # (throwing out all but the first at any single location),
                # so STR_{pos} is a unique ID
                finemap_input_z.write(
                    f'STR_{result_pos} {result_chrom:02} {result_pos} nan nan nan {beta} {se}\n'
                )
        if not any_strs:
            pathlib.Path(f"{workdir}/no_strs").touch()
            readme.write(
                "No nominally significant (p<=0.05) STRs were found in the region, "
                "so finemapping is being skipped.\n"
            )
            print(
                "No nominally significant (p<=0.05) STRs were found in the region, "
                "so finemapping is being skipped.",
                flush = True
            )
            return

        # load SNPs
        snps_to_filter = set()
        with open(filter_set_fname) as filter_file:
            next(filter_file) # skip header
            for line in filter_file:
                pos, ref, alt = line.strip().split('\t')[3:6]
                snps_to_filter.add((int(pos), ref, alt))

        with open(plink_results_fname) as plink_result_file:
            plink_results_reader = csv.reader(plink_result_file, delimiter='\t')
            header = next(plink_results_reader)
            cols = {
                col: header.index(col) for col in
                ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'BETA', 'SE', 'P', 'ERRCODE']
            }

            for result in plink_results_reader:
                result_chrom = int(result[cols['#CHROM']])
                result_pos = int(result[cols['POS']])
                ref = result[cols['REF']]
                alt = result[cols['ALT']]
                if (result_chrom, result_pos) < (chrom, start_pos):
                    continue
                if (result_chrom, result_pos) > (chrom, end_pos):
                    break
                if (result_pos, ref, alt) in snps_to_filter:
                    continue
                if result[cols['ERRCODE']] != '.' or float(result[cols['P']]) >= inclusion_threshold:
                    continue
                # snps can only be uniquely identified by pos, ref, alt
                # some IDs are duplicate at the same or even different locations
                # so those don't help

                beta = result[cols['BETA']]
                se = result[cols['SE']]
                finemap_input_z.write(
                    f'SNP_{result_pos}_{ref}_{alt} {result_chrom:02} {result_pos} {ref} {alt} nan {beta} {se}\n'
                )

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('phenotype')
    parser.add_argument('chrom', type=int)
    parser.add_argument('start_pos', type=int)
    parser.add_argument('end_pos', type=int)
    args = parser.parse_args()

    phenotype = args.phenotype
    chrom = args.chrom
    start_pos = args.start_pos
    end_pos = args.end_pos
    assert start_pos < end_pos

    outdir = f'{ukb}/finemapping/finemap_results/{phenotype}/{chrom}_{start_pos}_{end_pos}'

    with file_utils.temp_dir('finemap_write_input_variants', args) as tempdir:
        with open(f'{tempdir}/README.txt', 'w') as readme:
            write_input_variants(tempdir, readme, phenotype, chrom, start_pos, end_pos)
        file_utils.move_files(tempdir, outdir)

if __name__ == '__main__':
    main()

