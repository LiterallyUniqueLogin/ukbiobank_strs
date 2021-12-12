#!/usr/bin/env python3

import argparse
import json

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('outtable')
    parser.add_argument('outreadme')
    parser.add_argument('summary_tables', help='json dict of phenotype to table')
    args = parser.parse_args()

    with open(args.outreadme, 'w') as readme:
        readme.write(
            'Anything in the summary tables (i.e. association p-value <= 5e-8) that has FINEMAP pcausal '
            '>= 0.5 and is <= 10bp from the nearest exon boundary of an exon of a transcript that '
            'has support >=2\n'
        )

    first_table = True
    with open(args.outtable, 'w') as outfile:
        for phenotype, summary_table_fname in json.loads(args.summary_tables).items():
            with open(summary_table_fname) as table:
                if first_table:
                    header = next(table)
                    split = header.strip().split('\t')
                    nearby_exons_idx = split.index('nearby_exons')
                    pcausal_idx = split.index('pcausal')
                    outfile.write('phenotype\t')
                    outfile.write(header)
                else:
                    next_line = next(table)
                    if next_line != header:
                        print(next_line, header)
                        assert False

                for line in table:
                    split = line.strip().split('\t')
                    if split[pcausal_idx] != 'NA' and float(split[pcausal_idx]) < 0.5:
                        continue
                    exons = split[nearby_exons_idx]
                    if exons == 'none':
                        continue
                    close = False
                    for exon in exons.split(','):
                        dist = int(exon.split(':')[0])
                        if dist <= 10:
                            close = True
                            break
                    if not close:
                        continue
                    outfile.write(f'{phenotype}\t')
                    outfile.write(line)
                first_table = False

if __name__ == '__main__':
    main()
