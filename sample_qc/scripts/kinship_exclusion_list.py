import argparse

import list_from_sample_qc as sqc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('run_name')
    args = parser.parse_args()
    sqc.produce_list_from_sample_qc(args.run_name,
                                    'kinship_exclusion', False, 21)


if __name__ == '__main__':
    main()
