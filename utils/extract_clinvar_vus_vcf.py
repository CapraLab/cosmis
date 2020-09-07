#!/usr/bin/env python3

import csv
import json
import allel
from argparse import ArgumentParser


def parse_cmd():
    """
    Parse and return command-line arguments.

    Returns
    -------
    Command-line arguments.

    """
    parser = ArgumentParser()
    parser.add_argument('-o', '--output', dest='output', required=True,
                        type=str, help='Output file to store the results.')
    parser.add_argument('-v', '--vcf', dest='vcf', required=True, type=str,
                        help='Input file containing VEP annotations.')
    parser.add_argument('-c', '--csq', dest='csq', required=True, type=str,
                        help='Format line of the consequence annotations from'
                             'Ensembl VEP.')
    args = parser.parse_args()
    # do any necessary check on command-line arguments here
    return args


def main():
    """
    Extract variants from ClinVar that satisfy a pre-defined set of criteria.

    Returns
    -------

    """
    # parse command-line arguments
    args = parse_cmd()

    # parser CSQ header
    with open(args.csq, 'rt') as ipf:
        csq_header = [h.lower() for h in ipf.readline().strip().split('|')]
    print('Given CSQ header:')
    print(csq_header)

    # read the vcf file
    print('Now reading the CSQ records from the VCF file:', args.vcf)
    call_set = allel.read_vcf(input=args.vcf, fields='*')
    print('CSQ records read successfully.')
    
    # extract relevant fields from the vcf file
    clinsigs = call_set['variants/CLNSIG']
    clinstatus = call_set['variants/CLNREVSTAT']
    csq_annotations = call_set['variants/CSQ']

    # set up relevant fields and qualifiers
    vep_csq_fields = ['symbol', 'gene', 'feature', 'protein_position', 
        'amino_acids', 'ccds', 'ensp', 'swissprot', 'clin_sig', 'sift', 'polyphen']
    is_vus = 'Uncertain_significance'
    
    # store all qualified variants in a list
    qualified_variants = []

    # process each variant
    print('Now processing each variant ...')
    for x, y, z in zip(clinsigs, clinstatus, csq_annotations):
        # split VEP annotations into a list
        csqs = z.split('|')
        # retrieve the consequence of the current variant
        snp_type = csqs[csq_header.index('consequence')]
        # retrieve the canonical status of the current variant
        canonical = csqs[csq_header.index('canonical')]
        # consider only missense variants mapped to canonical transcripts
        if snp_type == 'missense_variant' and canonical == 'YES':
            qualified_variant = []
            # retrieved the CSQ fields that are of interest
            if x == is_vus:
                print(x, y, z)
                for k in vep_csq_fields:
                    qualified_variant.append(csqs[csq_header.index(k)])
                qualified_variants.append(qualified_variant) 
                print('Added', qualified_variant, 'to the list of qualified variants.')
        else:
            continue
                
    # write qualified variants to a csv file, one variant per row
    with open(args.output, 'wt') as csv_handle:
        csvwriter = csv.writer(csv_handle, delimiter=',')
        # write header row
        output_headers = vep_csq_fields + ['clinvar']
        csvwriter.writerow(output_headers)
        # write a row for each qualified variant
        for qualified_variant in qualified_variants:
            csvwriter.writerow(qualified_variant)
            

if __name__ == '__main__':
    main()