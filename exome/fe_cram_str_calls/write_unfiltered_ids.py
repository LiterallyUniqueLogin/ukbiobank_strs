import os

ukb = os.environ['UKB']

filtered_ids = set()

with open('filtered_ids.txt') as filtered_ids_file:
    for line in filtered_ids_file:
        line = line.strip()
        if '#' in line or line == '':
            continue
        filtered_ids.add(line)

ids = []

with open(f'{ukb}/side_analyses/exome_strs/snpstr_exome_str_cleaned_ids.txt') as ids_file:
    for line in ids_file:
        line = line.strip()
        if line == '' or line in filtered_ids:
            continue
        print(line)
