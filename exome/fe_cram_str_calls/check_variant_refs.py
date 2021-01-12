import gzip
import os

ukb = os.environ['UKB']

filtered_IDs = set()
with open('filtered_ids.txt') as filtered:
    for line in filtered:
        if '#' in line:
            continue
        line = line.strip()
        if line == '':
            continue
        filtered_IDs.add(line)

header = False
with open(f'{ukb}/side_analyses/exome_strs/snpstr_exome_str_cleaned_ids.txt') as ids:
    for num, ID in enumerate(ids):
        ID = ID.strip()
        if ID in filtered_IDs:
            continue

        reps = set()
        for batch in range(1, 11):
            with gzip.open(f'vcfs/{ID}_batch_{batch}.vcf.gz') as call_file:
                lines = call_file.readlines()
                reps.add(tuple(lines[-1].split()[:4]))
        if len(reps) != 1:
            if not header:
                print("Found STRs with different representations!")
                print("--------------------")
                header = True
            for rep in reps:
                print(rep)
            print("---------------------")
        print(f"Working on ID {num}", end="\r", flush=True)

if header:
    print("Done                     ")
else:
    print("All variants have some position and ref allele across batches")


