#!/usr/bin/env python3

import collections
import os

import cyvcf2

import str_utils
import python_utils as utils

ukb = os.environ['UKB']

refs = utils.load_reference(f'{ukb}/../../resources/dbase/human/hg19/hg19.fa')

str_loci = set()
with open(f'{ukb}/snpstr/str_loci.txt') as loci_file:
    for line in loci_file:
        str_loci.add(tuple(int(x) for x in line.split()))

with open(f'{ukb}/snpstr/repeat_units.tab', 'w') as out:
    out.write('chrom\tsnpstr_pos\tperiod\tunit\tcompound_unit\n')
    for chrom in range(1, 23):
        ref_chrom = refs[chrom].upper()
        vcf = cyvcf2.VCF(f'{ukb}/snpstr/info_field/chr{chrom}.vcf.gz')
        for var in vcf:
            if (int(var.CHROM), var.POS) not in str_loci:
                continue
            ref_seq = ref_chrom[(var.INFO['START'] - 1):var.INFO['END']]
            unit = str_utils.infer_repeat_unit(ref_seq, var.INFO['PERIOD'])
           
            compound_unit = False
            if unit is None:
                compound_unit = 'NA'
            else:
                if var.INFO['PERIOD'] % 2 == 0:
                    compound_unit = unit[:(len(unit)//2)] == unit[(len(unit)//2):]
                if var.INFO['PERIOD'] == 6:
                    compound_unit |= unit[:2] == unit[2:4] == unit[4:]
                
            out.write(f'{var.CHROM}\t{var.POS}\t{var.INFO["PERIOD"]}\t{unit}\t{compound_unit}\n')
