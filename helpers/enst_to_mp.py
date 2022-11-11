#!/usr/bin/env python3

import os, csv
import gzip
import json
import logging
from argparse import ArgumentParser
from Bio import SeqIO
from Bio import BiopythonWarning
import warnings
from cosmis import seq_utils

warnings.simplefilter('ignore', BiopythonWarning)


def parse_cmd():
    """
    Parse command-line arguments.

    Returns
    -------
    ArgumentParser
        Parsed command-line arguments.

    """
    parser = ArgumentParser()
    parser.add_argument('-c', '--config', dest='config', required=True,
                        type=str, help='A JSON file specifying options.')
    parser.add_argument('-t', '--transcripts', dest='transcripts', type=str,
                        required=True, help='''A list of ENSEMBL transcript IDs, 
                        one per line.''')
    parser.add_argument('-o', '--output', dest='output', required=False,
                        type=str, help='''Output file to which transcript 
                        mutation probabilities and mutation counts values will 
                        be written.''')
    parser.add_argument('-w', '--overwrite', dest='overwrite', required=False,
                        action='store_true', help='''Whether to overwrite already
                        computed MTR3D scores.''')
    parser.add_argument('-v', '--verbose', dest='verbose', required=False,
                        action='store_true', help='''Whether to output verbose 
                        data: number of contacting residues and number of 
                        missense and synonymous variants in the neighborhood of 
                        the  mutation site.''')
    return parser.parse_args()


def get_ensembl_accession(record):
    """

    Parameters
    ----------
    record : str
        Record ID is of the format: ">CCDS2.2|Hs109|chr1"

    Returns
    -------

    """
    parts = record.id.split('.')
    return parts[0]


def get_ccds_accession(record):
    """

    Parameters
    ----------
    record

    Returns
    -------

    """
    parts = record.id.split('|')
    return parts[0]


def parse_config(config):
    """

    Parameters
    ----------
    config

    Returns
    -------

    """
    with open(config, 'rt') as ipf:
        configs = json.load(ipf)

    # do necessary sanity checks before return
    return configs


def count_variants(variants):
    """
    Collects the statistics about position-specific counts of missense and
    synonymous variants.

    Parameters
    ----------
    variants : list
        A list of variant identifiers: ['A123B', 'C456D']

    Returns
    -------
    tuple
        Total number of synonymous and nonsynonymous variants.

    """
    missense_counts = 0
    synonymous_counts = 0
    for variant in variants:
        vv, _, _ = variant
        # skip variants whose MAF > 0.1%
        # if int(ac) / int(an) > 0.001:
        #    continue
        w = vv[0]  # wild-type amino acid
        v = vv[-1]  # mutant amino acid
        if w != v:  # missense variant
            missense_counts += 1
        else:  # synonymous variant
            synonymous_counts += 1
    return synonymous_counts, missense_counts


def main():
    """

    Returns
    -------

    """
    # configure the logging system
    logging.basicConfig(
        filename='enst_to_mp.log',
        level=logging.INFO,
        filemode='w',
        format='%(levelname)s:%(asctime)s:%(message)s'
    )

    # parse command-line arguments
    args = parse_cmd()

    # parse configuration file
    configs = parse_config(args.config)
    logging.info('Supplied configuration:')
    logging.info(json.dumps(configs, sort_keys=True, indent=4))

    # ENSEMBL cds
    print('Reading ENSEMBL CDS database ...')
    with gzip.open(configs['ensembl_cds'], 'rt') as cds_handle:
        ensembl_cds_dict = SeqIO.to_dict(
            SeqIO.parse(cds_handle, format='fasta'),
            key_function=get_ensembl_accession
        )

    # CCDS concensus coding sequences
    print('Reading NCBI CCDS database ...')
    with gzip.open(configs['ccds_cds'], 'rt') as ccds_handle:
        ccds_dict = SeqIO.to_dict(
            SeqIO.parse(ccds_handle, format='fasta'),
            key_function=get_ccds_accession
        )

    # parse gnomad transcript-level variants
    print('Reading gnomAD variant database ...')
    with open(configs['gnomad_variants'], 'rt') as variant_handle:
        # transcript_variants will be a dict of dicts where major version
        # ENSEMBL transcript IDs are the first level keys and "ccds", "ensp",
        # "swissprot", "variants" are the second level keys. The value of each
        # second-level key is a Python list.
        transcript_variants = json.load(variant_handle)

    output_dir = os.path.abspath(configs['output_dir'])

    # compute the mutation probabilities and variant counts for each transcript
    mutation_prob_vs_count = []
    with open(args.transcripts, 'rt') as ipf:
        for transcript in ipf:
            transcript = transcript.strip()
            print('Processing transcript %s' % transcript)

            # get the coding sequence of the transcript
            try:
                transcript_cds = ensembl_cds_dict[transcript].seq
            except KeyError:
                print('No CDS found in Ensembl CDS database! Looking for it '
                      'in the CCDS database ...')
                transcript_cds = None

            if transcript_cds is None:
                try:
                    ccds_id = transcript_variants[transcript]['ccds'][0]
                    transcript_cds = ccds_dict[ccds_id].seq
                except KeyError:
                    print('No CDS found in CCDS database! Skipped.')
                    continue

            # check that the CDS does not contain invalid nucleotides
            if not seq_utils.is_valid_cds(transcript_cds):
                print('Invalid CDS! Skipped.')
                continue

            # determine the total rate for the current transcript
            syn_prob, mis_prob = seq_utils.get_transcript_mutation_prob(
                transcript_cds
            )

            # get all variants of this transcript reported in gnomAD
            try:
                variants = transcript_variants[transcript]['variants']
            except KeyError:
                logging.critical('No variants found in %s in gnomAD',
                                 transcript)
                logging.critical('%s was skipped ...', transcript)
                continue

            # tabulate variants at each site
            syn_count, mis_count = count_variants(variants)
            mutation_prob_vs_count.append([
                transcript,
                len(transcript_cds) // 3 - 1,  # -1 to subtract stop codon
                '{:.3e}'.format(syn_prob),
                syn_count,
                '{:.3e}'.format(mis_prob),
                mis_count
            ])

        with open(file=os.path.join(output_dir, args.output), mode='wt') as opf:
            header = [
                'enst_id', 'length', 'syn_prob', 'syn_count', 'mis_prob', 'mis_count'
            ]
            csv_writer = csv.writer(opf, delimiter='\t')
            csv_writer.writerow(header)
            csv_writer.writerows(mutation_prob_vs_count)


if __name__ == '__main__':
    main()
