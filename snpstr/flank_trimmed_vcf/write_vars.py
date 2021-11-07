#!/usr/bin/env python3

import cyvcf2

with open('vars.tab', 'w') as out:
    out.write('chrom\tpos\tend_pos\tsnpstr_pos\n')
    for chrom in range(1, 23):
        vcf = cyvcf2.VCF(f'chr{chrom}.vcf.gz')
        for var in vcf:
            out.write(
                f'{chrom}\t{var.POS}\t{var.POS + len(var.REF) - 1}\t{var.INFO["ORIGINAL_POS"]}\n'
            )
