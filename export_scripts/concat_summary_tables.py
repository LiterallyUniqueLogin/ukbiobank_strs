#!/usr/bin/env python3

import argparse
import json

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('outfile')
    parser.add_argument('summary_tables', help='json dict of phenotype to table')
    args = parser.parse_args()

    with open(args.outfile, 'w') as outfile:
        first_table = True
        for phenotype, summary_table_fname in json.loads(args.summary_tables).items():
            with open(summary_table_fname) as table:
                if first_table:
                    header = next(table)
                    outfile.write('phenotype\t')
                    outfile.write(header)
                else:
                    next_line = next(table)
                    if next_line != header:
                        print(next_line, header)
                        assert False

                for line in table:
                    outfile.write(f'{phenotype}\t')
                    outfile.write(line)
            first_table = False

if __name__ == '__main__':
    main()
