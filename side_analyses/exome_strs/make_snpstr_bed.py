import gzip
import os
import re
import sys

ukb = os.environ['UKB']

str_ids = set()
with open(f'{ukb}/snpstr/str_ids.txt') as str_ids_file:
    for line in str_ids_file:
        str_ids.add(line.strip())

# based on HipSTR source code it looks like number of copies field in the
# bed is ignored
# https://github.com/tfwillems/HipSTR/blob/master/src/region.cpp
with open('snpstr_strs_19.bed', 'w') as bed:
    start_re = re.compile('START=([0-9]+)')
    end_re = re.compile('END=([0-9]+)')
    period_re = re.compile('PERIOD=([0-9]+)')
    for chrom in range(1, 23):
        with gzip.open(f'{ukb}/snpstr/info_field/chr{chrom}.vcf.gz', 'rt') as info_vcf:
            for line in info_vcf:
                if line[0] == '#':
                    continue
                split = line.split()
                chrom_str, _, name, ref = split[:4]
                if name not in str_ids:
                    continue
                str_ids.remove(split[2])
                start = int(start_re.search(line).group(1)) - 1
                end = end_re.search(line).group(1)
                period = period_re.search(line).group(1)
                bed.write(f'chr{chrom_str}\t{start}\t{end}\t{period}\t{len(ref)/int(period)}\t{name}\n')
        print(f'Done with chr{chrom}', flush=True, file=sys.stderr)

assert len(str_ids) == 0
