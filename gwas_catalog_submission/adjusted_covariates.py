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
        curr_out = 'age|sex|genetic PCs 1-40|'
        assessments = line[3].split(',')
        if len(assessments) == 1:
            assert assessments[0] == 'initial_assessment'
            curr_out += 'Excluding participants without measurements from the initial assessment, no covariate needed'
        else:
            assert len(assessments) == 2
            curr_out += 'Using a categorical covariate to distinguish between participants with initial and repeat assessment measurements'
            if 'device_id' in line[4]:
                curr_out += ', excluding other participants'
        if line[4] != 'None':
            if ',' not in line[4]:
                curr_out += f'|Using categorical covariates to distinguish between {line[4]}s'
            else:
                covar, val_str = line[4].replace(' ', '')[:-1].split('(')
                covar = covar.replace(',', '')
                vals = val_str.split(',')
                curr_out += '|Using '
                if len(vals) == 2:
                    curr_out += 'a categorical covariate'
                else:
                    curr_out += 'categorical covariates'
                curr_out += f' to distinguish between {covar}s {val_str}, excluding participants with other {covar}s'

        output[(line[0], line[1])] = curr_out

with open('supp_table_1.csv') as wb_file:
    first = True
    for line in csv.reader(wb_file):
        if first:
            first = False
            continue
        curr_out = 'age|sex|genetic PCs 1-40|Using a categorical covariate to distinguish between participants with initial and repeat assessment measurements'
        if line[8] != 'NA':
            curr_out += ', excluding other participants'
        covar_tuple = ast.literal_eval(line[12])
        if covar_tuple is not None:
            if len(covar_tuple) == 2:
                assert 'device_id' in covar_tuple[0]
                curr_out += '|Using categorical covariates to distinguish between device_ids'
            else:
                assert 'aliquot' in covar_tuple[0]
                curr_out += f'|Using categorical covariates to distinguish between aliquots {",".join(str(num) for num in covar_tuple[2])}, excluding participants with other aliquots'

        output[('white_british', line[0])] = curr_out

for key in sorted(output.keys()):
    print(output[key])
