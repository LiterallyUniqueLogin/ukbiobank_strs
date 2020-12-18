import argparse

parser = argparse.ArgumentParser()
parser.add_argument('summary_file')
args = parser.parse_args()

filtered_IDs = set()
with open('filtered_ids.txt') as filtered:
    for line in filtered:
        if '#' in line:
            continue
        line = line.strip()
        if line == '':
            continue
        filtered_IDs.add(line)


with open(args.summary_file) as summary:
    prev_IDs = {}
    for line in summary:
        if '_batch_' in line:
            ID, _ = line.split('_batch_')
            if ID in filtered_IDs:
                continue
            if ID in prev_IDs:
                prev_IDs[ID] += 1
            else:
                prev_IDs[ID] = 1
        else:
            if len(prev_IDs) > 0:
                for ID, count in prev_IDs.items():
                    print(f'({count}) {ID}')
            prev_IDs = {}
            print(line, end='')

