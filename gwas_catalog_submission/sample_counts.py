#!/usr/bin/env python3

import ast
import csv

output = {}
with open('supp_table_9.csv') as eth_file:
    first = True
    for line in csv.reader(eth_file, delimiter=',', quotechar='"'):
        if first:
            first = False
            continue
        output[(line[0], line[1])] = line[2]

with open('supp_table_1.csv') as wb_file:
    first = True
    for line in csv.reader(wb_file):
        if first:
            first = False
            continue
        output[('white_british', line[0])] = line[11]

for key in sorted(output.keys()):
    print(output[key])
