#!/usr/bin/env python3

"""
@TODO A description of this script needs to be added.
"""

import json
# from Bio import bgzf
import gzip
from collections import defaultdict
from argparse import ArgumentParser


def parse_cmd():
    """
    Parses command-line arguments.

    Returns
    -------
    ArgumentParser
        An ArgumentParser object containing parsed command-line arguments.

    """
    parser = ArgumentParser()
    parser.add_argument(
        '-o', '--output', dest='output', required=True, type=str,
        help='Output file to store the results.'
    )
    parser.add_argument(
        '-i', '--input', dest='vcf', required=True, type=str,
        help='Input VCF file.'
    )
    parser.add_argument(
        '-v', '--vep-header', dest='vep', required=True, type=str,
        help='Format line of the consequence annotations from Ensembl VEP.'
    )
    args = parser.parse_args()
    # do any necessary check on command-line arguments here
    return args


def main():
    # parse command-line arguments
    args = parse_cmd()

    # parse vep header
    with open(args.vep, 'rt') as ipf:
        vep_header = [h.lower() for h in ipf.readline().strip().split('|')]
    print('Given VEP header:')
    print(vep_header)

    #
    print('Now reading the VEP records from the VCF file:', args.vcf)
    vcf_header = [
        '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'
    ]
    transcript_variants = defaultdict(lambda: defaultdict(list))
    with gzip.open(args.vcf, 'rt') as ipf:
        for l in ipf:
            # skip comment lines
            if l.startswith('#'):
                continue
            variant_fields = dict(zip(vcf_header, l.strip().split('\t')))

            # skip variant calls that does not meet quality requirement
            if variant_fields['FILTER'] != 'PASS':
                continue

            # get the INFO field
            info_fields = dict(
                item.split('=') for item in variant_fields['INFO'].split(';')
                if '=' in item
            )

            # get VEP annotations
            vep_annotations = info_fields['vep'].split(',')

            # VEP records that passed the filters
            for annot in vep_annotations:
                annot = annot.split('|')
                if annot[vep_header.index('biotype')] != 'protein_coding':
                    # print('Skipped a record of non-protein-coding variant.')
                    continue
                consequence = annot[vep_header.index('consequence')]
                transcript_id = annot[vep_header.index('feature')]
                ccds_id = annot[vep_header.index('ccds')]
                ensp_id = annot[vep_header.index('ensp')]
                swissprot_id = annot[vep_header.index('swissprot')]
                protein_position = annot[vep_header.index('protein_position')]
                amino_acids = annot[vep_header.index('amino_acids')]

                # add current record
                transcript_variants[transcript_id]['ccds'].append(ccds_id)
                transcript_variants[transcript_id]['ensp'].append(ensp_id)
                transcript_variants[transcript_id]['swissprot'].append(swissprot_id)

                #
                if consequence == 'missense_variant':
                    w, v = amino_acids.split('/')
                    variant = w + protein_position + v
                elif consequence == 'synonymous_variant':
                    variant = amino_acids + protein_position + amino_acids
                else:
                    continue

                transcript_variants[transcript_id]['variants'].append(
                    (variant, int(info_fields['AC']), int(info_fields['AN']))
                )

    # remove duplicate ids
    for k in transcript_variants.keys():
        transcript_variants[k]['ccds'] = list(set(transcript_variants[k]['ccds']))
        transcript_variants[k]['ensp'] = list(set(transcript_variants[k]['ensp']))
        transcript_variants[k]['swissprot'] = list(set(transcript_variants[k]['swissprot']))

    # serialize variants into JSON format
    print('Now serializing variants in JSON format ...')
    with open(args.output, 'wt') as json_handle:
        json.dump(transcript_variants, json_handle, sort_keys=True, indent=4)


if __name__ == '__main__':
    main()
