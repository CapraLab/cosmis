#!/usr/bin/env python3

import json
from argparse import ArgumentParser
from collections import defaultdict


def parse_cmd_args():
    """
    Parse command-line arguments.

    Returns
    -------
    Parsed command-line arguments.

    """
    parser = ArgumentParser(description='')
    parser.add_argument(
        '-i', '--input', dest='input', type=str, required=True,
        help='''Input file containing gnomAD sequencing depths of coverage.'''
    )
    parser.add_argument(
        '-o', '--output', dest='output', type=str, required=False,
        default='out.json', help='''Output file.'''
    )
    # do any necessary checks on command-line arguments here
    return parser.parse_args()


def main():
    """

    Returns
    -------

    """
    # parse all command-line arguments
    cmd_args = parse_cmd_args()

    #
    coord_to_seqcov = defaultdict(dict)
    with open(cmd_args.input, 'rt') as ip_handle:
        # remove the header line
        _ = next(ip_handle)
        # process each position
        prev_chrom = 'chr1'
        print('Now processing chr1 ...')
        for line in ip_handle:
            # chrom, start, end, mean_cov, median_cov, ...
            chrom, _, pos, mean, median = line.strip().split()
            if chrom != prev_chrom:
                print('Now processing {} ...'.format(chrom))
                prev_chrom = chrom
            coord_to_seqcov[chrom][int(pos)] = (float(mean), float(median))


    print('Now writing depths of coverage to {}'.format(cmd_args.output))
    with open(cmd_args.output, 'wt') as op_handle:
        json.dump(coord_to_seqcov, fp=op_handle, indent=4)


if __name__ == '__main__':
    main()
