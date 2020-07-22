"""
Extract a subset of fields from the ukb29170.csv file by name.

Always also extracts the id field.
Output form is csv
"""

import csv
import os
import sys

if len(sys.argv) < 2:
    print("Excpecting at least one argument: the column header(s) to look for")
    sys.exit(-1)

fields = sys.argv[1:]
fields.insert(0, "eid")

ukb = os.environ['UKB']

def print_line(idxs, line):
    for count, idx in enumerate(idxs):
        if count != len(idxs) - 1:
            print(line[idx], end=",")
        else:
            print(line[idx])


with open(f'{ukb}/main_dataset/raw_data/ukb29170.csv') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar="\"")
    header = next(reader)
    field_idxs = []
    for field in fields:
        # this will raise an error if the field isn't found
        field_idxs.append(header.index(field))

    # print the header
    print_line(field_idxs, header)

    # Extract the info
    for line in reader:
        print_line(field_idxs, line)

