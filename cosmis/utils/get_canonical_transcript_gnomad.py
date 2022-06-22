#!/usr/bin/env python3

import gzip
from argparse import ArgumentParser


def parse_cmd():
    """

    Returns
    -------

    """
    parser = ArgumentParser()
    parser.add_argument('-o', '--output', dest='output', required=True,
                        type=str, help='Output file to store the results.')
    parser.add_argument('-i', '--input', dest='vcf', required=True, type=str,
                        help='Input VCF file.')
    parser.add_argument('-v', '--vep-header', dest='vep', required=True, type=str,
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

    # parser CSQ header
    with open(args.vep, 'rt') as ipf:
        vep_header = [h.lower() for h in ipf.readline().strip().split('|')]
    print('Given VEP header:')
    print(vep_header)

    # index of canonical field
    canonical_index = vep_header.index('canonical')
    transcript_index = vep_header.index('feature')
    ensg_index = vep_header.index('gene')

    #
    print('Now reading the VEP records from the VCF file:', args.vcf)
    # Biopython's bgzf module for reading bgz files
    vcf_header = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

    # canonical transcript
    canonical_transcripts = set()

    # with bgzf.open(args.vcf, 'rt') as ipf:
    with gzip.open(args.vcf, 'rt') as ipf:
        for l in ipf:
            if l.startswith('#'):
                continue
            variant_fields = dict(zip(vcf_header, l.strip().split('\t')))

            # skip variant calls that does not meet quality requirement
            if variant_fields['FILTER'] != 'PASS':
                continue

            # get the INFO field
            info_fields = dict(item.split('=') for item in
                               variant_fields['INFO'].split(';') if '=' in item)

            # get VEP annotations
            vep_annotations = info_fields['vep'].split(',')

            # CSQs that passed the filters
            for annot in vep_annotations:
                annot = annot.split('|')
                if annot[vep_header.index('biotype')] != 'protein_coding':
                    continue
                if annot[canonical_index] == 'YES':
                    canonical_transcripts.add((annot[ensg_index], annot[transcript_index]))

    # serialize variants into JSON format
    print('Now writing canonical transcripts to file ...')
    with open(args.output, 'wt') as opf:
        for pair in canonical_transcripts:
            opf.write('\t'.join(pair))
            opf.write('\n')


if __name__ == '__main__':
    main()
