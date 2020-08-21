"""
Exome files include all STR variants in the SNPSTR vcfs that intersect with the
exonic regions even if the intersection is only flanking indels
Instead, filter to STRs that are wholly contained within the exonic region

cat chr*tab | awk '{print $1 " " $2}' | sort | uniq -c | sort
showed that there was a duplication variant 2:96780975, I manually deleted it
before running this script

Counting the number of lines in the output in this directory compared to the
output in $UKB/side_analyses/entropy/exome_filtered shows that there are the
same number of lines (ignore the 22 header lines there), these ones are just
the canonical locations without the nearby SNPs being included in the variant
"""

import os
import subprocess as sp
import sys

ukb = os.environ['UKB']

def main():
    filedir = f'{ukb}/misc_data/hipstr/bed'

    for chrom in range(1, 23):
        chrom = str(chrom)
        with open(f'{filedir}/hipstr_snpstr_exome_chr{chrom}.bed') as infile, \
             open(f'{filedir}/filtered_hipstr_snpstr_exome_chr{chrom}.bed', 'w') as outfile:
            nvariants = 0
            for inline in infile:
                nvariants += 1
                print(f"Working on chr {chrom} variant {nvariants}",
                      end='\r')
                str_start, str_end = inline.split()[1:3]
                str_start = int(str_start)
                str_end = int(str_end)
                with open(f'{ukb}/side_analyses/exome_strs/hg19_final.bed') as exome_bed:
                    for rline in exome_bed:
                        rchrom, rstart, rend = rline.split()
                        rstart = int(rstart)
                        rend = int(rend)
                        if rchrom != chrom:
                            continue
                        if rend > str_start:
                            if rstart < str_start and str_end < rend:
                                outfile.write(inline)
                            break

if __name__ == "__main__":
    main()
