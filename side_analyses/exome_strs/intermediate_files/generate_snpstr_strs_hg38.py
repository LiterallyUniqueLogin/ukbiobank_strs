import os

ukb = os.environ['UKB']

str_to_line = {}
with open(f'{ukb}/side_analyses/exome_strs/hg38.hipstr_reference.bed') as ref:
    for linenum, line in enumerate(ref):
        str_to_line[line.split()[5]] = line
        if linenum % 1000 == 0:
            print(f'Loaded {linenum} lines', end='\r', flush=True)

already_seen_ids = set()
with open(f'{ukb}/snpstr/str_ids.txt') as str_ids, \
        open(f'{ukb}/side_analyses/exome_strs/snpstr_strs_hg38.bed', 'w') as snpstr_strs_hg38:
    for line in str_ids:
        str_id = line.strip()
        if str_id in already_seen_ids:
            continue
        already_seen_ids.add(str_id)
        if f'Human_{str_id}' in str_to_line:
            snpstr_strs_hg38.write(str_to_line[f'Human_{str_id}'])
        else:
            snpstr_strs_hg38.write(str_to_line[str_id])


