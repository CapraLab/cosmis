#!/usr/bin/env python3

import json
from argparse import ArgumentParser


def parse_cmd():
    """

    Returns
    -------

    """
    parser = ArgumentParser()
    parser.add_argument('-f', '--filters', dest='filters', required=True,
                        type=str, help='JSON file specifying the filters')
    parser.add_argument('-o', '--output', dest='output', required=True,
                        type=str, help='Output file to store the results.')
    parser.add_argument('-v', '--vcf', dest='vcf', required=True, type=str,
                        help='Input VCF file.')
    parser.add_argument('-c', '--csq', dest='csq', required=True, type=str,
                        help='Format line of the consequence annotations from'
                             'Ensembl VEP.')
    args = parser.parse_args()
    # do any necessary check on command-line arguments here
    return args


def main():
    """

    Returns
    -------

    """
    # parse command-line arguments
    args = parse_cmd()

    # parse filters
    with open(args.filters, 'rt') as json_file:
        filters = json.load(json_file)
    print('Given filters are:')
    print(filters)

    # parser CSQ header
    with open(args.csq, 'rt') as ipf:
        csq_header = [h.lower() for h in ipf.readline().strip().split('|')]
    print('Given CSQ header:')
    print(csq_header)

    #
    print('Now reading the CSQ records from the VCF file:', args.vcf)
    call_set = allel.read_vcf(input=args.vcf, fields='variants/CSQ')
    csq_annotations = call_set['variants/CSQ']
    print('CSQ records read successfully.')

    # CSQs that passed the filters
    print('Now filtering the CSQ records')
    with open(args.output, 'wt') as opf:
        opf.write(','.join(csq_header) + '\n')
        for csq in csq_annotations:
            csq = csq.split('|')
            passed = True
            for k, v in filters.items():
                if csq[csq_header.index(k)] != v:
                    print('Skip current CSQ record because the', k, 'filter is not met.')
                    passed = False
                    break
            if passed:
                opf.write(','.join(csq) + '\n')


if __name__ == '__main__':
    main()
