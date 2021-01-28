import os
import subprocess as sp
import time

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
        ids.append(line)

total_time = 0
for id_count, curr_id in enumerate(ids):
    id_count += 1
    print(f"Working on id #{id_count}   ", end='', flush=True)
    files = []
    for i in range(1, 11):
        files.append(f'{ukb}/exome/fe_cram_str_calls/vcfs/{curr_id}_batch_{i}.vcf.gz')
    files = ','.join(files)

    command = f"""
    source ~/.bashrc ;
    conda activate ukb ;
    cd {ukb}/../trtools/repo ;
    python -m trtools.mergeSTR.mergeSTR \
            --vcfs {files} \
            --vcftype hipstr \
            --out {ukb}/exome/fe_cram_str_calls/merged_vcfs/{curr_id}.merged \
            --trim
    """

    start = time.time()
    out = sp.run(command, shell=True, capture_output=True)
    total_time += time.time() - start
    avg_time = total_time/id_count
    ttd = avg_time*(630 - len(filtered_ids) - id_count)
    print(f'{avg_time:0.1f}sec per TR. Time till done: {ttd:0.1f}', end='\r',
         flush=True)
    if out.returncode != 0:
        print(f'STR: {curr_id}')
        print('-----stdout-----')
        print(out.stdout.decode())
        print('-----stderr-----')
        print(out.stderr.decode())
        out.check_returncode()
