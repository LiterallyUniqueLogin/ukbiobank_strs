"""
Exome files include all STR variants in the SNPSTR vcfs that intersect with the
exonic regions even if the intersection is only flanking indels
Instead, filter to STRs that are wholly contained within the exonic region

cat chr*tab | awk '{print $1 " " $2}' | sort | uniq -c | sort
demonstrates that all exome starting positions are correctly unique
so I can index off of that
"""

import argparse
import os
import subprocess as sp
import sys

ukb = os.environ['UKB']

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--het", action="store_true")

    args = parser.parse_args()

    indir = f'{ukb}/side_analyses/entropy/exome'
    outdir = f'{ukb}/side_analyses/entropy/exome_filtered'
    if args.het:
        indir += '_het'
        outdir += '_het'

    sp.run(f"cp {indir}/README.txt {outdir}", shell=True)

    for chrom in range(1, 23):
        chrom = str(chrom)
        with open(f'{indir}/chr{chrom}.tab') as infile, \
             open(f'{outdir}/chr{chrom}.tab', 'w') as outfile:
            first = True
            nvariants = 0
            for inline in infile:
                if first:
                    # header
                    outfile.write(inline)
                    first = False
                    continue
                nvariants += 1
                print(f"Working on chr {chrom} variant {nvariants}",
                      end='\r')
                pos = inline.split()[1]
                command = f'''
                    source ~/.bashrc ;
                    conda activate ukb ; 
                    bcftools query -f '%INFO/START %INFO/END\\n' -r {chrom}:{pos} \\
                        {ukb}/snpstr/info_field/chr{chrom}.vcf.gz 
                    '''
                out = sp.run(command, shell=True, capture_output=True, text=True)
                if out.stderr:
                    print(out.stderr)
                    sys.exit()
                str_start, str_end = out.stdout.split()
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
