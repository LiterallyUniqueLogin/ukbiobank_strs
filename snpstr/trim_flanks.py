#!/usr/bin/env python3

import gzip
import os
import subprocess as sp

ukb = os.environ['UKB']

def trim_flank(chrom):
    with gzip.open(f'{ukb}/snpstr/info_field/chr{chrom}.vcf.gz', 'rt') as info_vcf, \
            open(f'{ukb}/snpstr/flank_trimmed_vcf/chr{chrom}.vcf', 'w') as trimmed_vcf:
        emitted_info = False
        for line in info_vcf:
            if line[0] == '#':
                if line[1] == '#':
                    if line[2:6] == 'INFO':
                        if not emitted_info:
                            trimmed_vcf.write("##INFO=<ID=ORIGINAL_POS,Number=1,Type=Integer,Description=\"The start position before flank trimming\">\n")
                            emitted_info = True
                        continue
                    else:
                        trimmed_vcf.write(line)
                else:
                    trimmed_vcf.write("##command=trim_flanks.py\n")
                    trimmed_vcf.write(line) # header
            else:
                split = line.strip().split()
                info_split = split[-1].split(';')
                start = None
                end = None
                for part in info_split:
                    if part.startswith('START='):
                        start = int(part[6:])
                    if part.startswith('END='):
                        end = int(part[4:])
                if end is None or start is None:
                    print(line)
                    print(info_split)
                pos = int(split[1])
                start_diff = start - pos
                new_end = end - pos + 1
                trimmed_vcf.write('\t'.join([
                    split[0],
                    str(start),
                    split[2],
                    split[3][start_diff:new_end],
                    '.',
                    split[5],
                    split[6],
                    f'ORIGINAL_POS={pos}'
                ]))
                trimmed_vcf.write('\n')

for chrom in range(1, 23):
    trim_flank(chrom)
    sp.run(
        f'bgzip -f {ukb}/snpstr/flank_trimmed_vcf/chr{chrom}.vcf && '
        f'tabix {ukb}/snpstr/flank_trimmed_vcf/chr{chrom}.vcf.gz',
        shell=True,
        check=True
    )


