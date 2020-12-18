import gzip
import sys

# confirm that the STRs named in the info vcfs correspond to the STRs named in
# the SNPSTR vcfs

str_ids = set()
with open("str_ids.txt") as str_ids_file:
    for line in str_ids_file:
        str_ids.add(line.strip())

missing_names = set()
found_str_ids = set()
prev_str_if_failed = None
for chrom in range(1, 23):
    # id -> pos, ref, alt
    info_str_by_id = {}
    with gzip.open(f"info_field/chr{chrom}.vcf.gz", 'rt') as info_vcf:
        for line in info_vcf:
            if line[0] == '#':
                continue
            split = line.split()
            info_str_by_id[split[2]] = (split[1], split[3], split[4])
    with open(f'vcf_1_sample/chr{chrom}.vcf') as snpstr_vcf:
        for line in snpstr_vcf:
            if line[0] == '#':
                continue
            split = line.split()
            name = split[2]
            if prev_str_if_failed is not None:
                if name != prev_str_if_failed[0]:
                    print(f'{prev_str_if_failed} in SNPSTR failed to match '
                          f'INFO vcf {info_str_by_id[prev_str_if_failed[0]]}')
                else:
                    print(f'Found a duplicate STR ID {name}')

            if name not in str_ids:
                prev_str_if_failed = None
                continue
            found_str_ids.add(name)
            split = line.split()
            snpstr_data = (split[1], split[3], split[4])
            if name not in info_str_by_id:
                missing_names.add(name)
                continue
            if info_str_by_id[name] != snpstr_data:
                prev_str_if_failed = (name, snpstr_data)
            else:
                prev_str_if_failed = None
    print(f'Done with chr {chrom}', file=sys.stderr)


print('------------')
if len(missing_names) == 0:
    print('All SNPSTR STR ids present in the INFO vcfs')
else:
    print('Some SNPSTR STR ids not present in INFO vcfs:', missing_names)

print('-------------')
if found_str_ids == str_ids:
    print('Found all str_ids')
else:
    print('str_ids not found', str_ids - found_str_ids)

