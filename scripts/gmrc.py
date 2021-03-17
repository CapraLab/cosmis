#!/usr/bin/env python3


##################################################
# This script computes gene-level mutation rates.
# Hence, the name gmrc.py
##################################################


import sys, os, re, csv
import gzip
import json
import logging
import urllib
import numpy as np
from collections import defaultdict
from argparse import ArgumentParser
from mtr3d.struct.contact import Contact, BACKBONE_ATOMS
from mtr3d.mapping.sifts import SIFTS
from mtr3d.mapping.ensembl_uniprot_pdb import EnsemblUniProtPDB
from mtr3d.utils import pdb_utils
from mtr3d.utils.genetic_code import GENETIC_CODE
from Bio import SeqIO
from Bio.SeqUtils import seq1
from Bio.PDB import PDBParser, NeighborSearch, is_aa, PDBList
from mtr3d.mutation_rates.trinucleotide_context_rates import MUTATION_RATES_UNIQUE

from Bio import BiopythonWarning
import warnings
warnings.simplefilter('ignore', BiopythonWarning)


def get_codon_mutation_rates(cds):
    """

    Parameters
    ----------
    cds : str
        Coding sequence.

    Returns
    -------
    mutation probabilities : list
        A list of tuples consisting of synonymous and nonsynonymous mutation 
        probabilities.
    """
    if len(cds) % 3 != 0:
        raise ValueError('Given CDS length is not a multiple of 3.')

    num_codons = len(cds) // 3

    mutation_rates = []
    # one nucleotide before and one nucleotide after the codon
    for codon_number in range(1, num_codons + 1):
        # determine the codon sequence
        codon_sequence = cds[(codon_number - 1) * 3:codon_number * 3]
        
        synonymous_rate = 0
        nonsynonymous_rate = 0
        
        # determine the mutation rate of the first and the last codons
        # consider only two mutatable nucleotides
        if codon_number == 1 or codon_number == num_codons:
            # first codon
            if codon_number == 1:
                sequence_context = cds[:4]
                # i is the zero-indexed position of the mutated nucleotide
                for i in range(1, 3):
                    trinucleotide = sequence_context[i - 1:i + 2]
                    rates = MUTATION_RATES_UNIQUE[trinucleotide]
                    for k, v in rates.items():
                        mutant_sequence = codon_sequence[:i] + k[1] + codon_sequence[i + 1:]
                        if GENETIC_CODE[codon_sequence] == GENETIC_CODE[mutant_sequence]:
                            synonymous_rate += v
                        elif GENETIC_CODE[codon_sequence] != GENETIC_CODE[mutant_sequence] \
                            and GENETIC_CODE[mutant_sequence] != 'STOP':
                            nonsynonymous_rate += v
            # last codon        
            else:
                sequence_context = cds[-4:]
                # i is the zero-indexed position of the mutated nucleotide
                for i in range(0, 2):
                    trinucleotide = sequence_context[i:i + 3]
                    rates = MUTATION_RATES_UNIQUE[trinucleotide]
                    for k, v in rates.items():
                        mutant_sequence = codon_sequence[:i] + k[1] + codon_sequence[i + 1:]
                        if GENETIC_CODE[codon_sequence] == GENETIC_CODE[mutant_sequence]:
                            synonymous_rate += v
                        elif GENETIC_CODE[codon_sequence] != GENETIC_CODE[mutant_sequence] \
                            and GENETIC_CODE[mutant_sequence] != 'STOP':
                            nonsynonymous_rate += v

        # codons other than the first and the last
        # consider all three mutatable nucleotides
        else:
            # one nucleotide before and one nucleotide after the codon
            sequence_context = cds[(codon_number - 1) * 3 - 1:(codon_number - 1) * 3 + 4]
            # mutate nucleotide in the codon iteratively
            for i in range(3):
                trinucleotide = sequence_context[i:i + 3]
                rates = MUTATION_RATES_UNIQUE[trinucleotide]
                for k, v in rates.items():
                    codon_seq_list = list(codon_sequence)
                    codon_seq_list[i] = k[1]
                    mutant_sequence = ''.join(codon_seq_list)
                    if GENETIC_CODE[codon_sequence] == GENETIC_CODE[mutant_sequence]:
                        synonymous_rate += v
                    elif GENETIC_CODE[codon_sequence] != GENETIC_CODE[mutant_sequence] \
                        and GENETIC_CODE[mutant_sequence] != 'STOP':
                        nonsynonymous_rate += v

        mutation_rates.append((synonymous_rate, nonsynonymous_rate))
    return mutation_rates


def parse_cmd():
    """

    Returns
    -------

    """
    parser = ArgumentParser()
    # parser.add_argument('-a', '--aa-seq', dest='aa_sequence', type=str, required=True,
    #                     help='FASTA file representing the UniProt canonical '
    #                          'amino acid sequence of the protein')
    # parser.add_argument('-n', '--non', dest='nonsynonymous', type=str, required=True,
    #                     help='A tab-delimited file in which each row specifies a'
    #                          'position-specific count of nonsynonymous variants')
    # parser.add_argument('-s', '--syn', dest='synonymous', type=str, required=True,
    #                     help='A tab-delimited file in which each row specifies a '
    #                          'position-specific count of synonymous variants.')
    parser.add_argument('-c', '--config', dest='config', required=True,
                        type=str, help='A JSON file specifying options.')
    parser.add_argument('-t', '--transcripts', dest='transcripts', type=str,
                        required=True, help='A list of ENSEMBL transcript IDs, '
                        'one per line.')
    parser.add_argument('-o', '--output', dest='output', required=False,
                        type=str, help='Output file to which site-specific MTR '
                        'values will be written.')
    parser.add_argument('-w', '--overwrite', dest='overwrite', required=False,
                        action='store_true', help='Whether to overwrite already '
                        'computed MTR3D scores.')
    parser.add_argument('-v', '--verbose', dest='verbose', required=False,
                        action='store_true', help='Whether to output verbose '
                        'data: number of contacting residues and number of '
                        'missense and synonymous variants in the neighborhood'
                        'of the mutation site.')
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
    #
    missense_counts = 0
    synonymous_counts = 0
    for variant in variants:
        vv, ac, an = variant
        # skip variants whose MAF > 0.01%
        if int(ac) / int(an) > 0.001:
            continue
        w = vv[0]  # wild-type amino acid
        v = vv[-1]  # mutant amino acid
        pos = vv[1:-1]  # position in the protein sequence
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
        filename='mtr3d.log',
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
        ensembl_cds_dict = SeqIO.to_dict(SeqIO.parse(cds_handle, format='fasta'),
                                 key_function=get_ensembl_accession)

    # CCDS concensus coding sequences
    print('Reading NCBI CCDS database ...')
    with gzip.open(configs['ccds_cds'], 'rt') as ccds_handle:
        ccds_dict = SeqIO.to_dict(SeqIO.parse(ccds_handle, format='fasta'),
                                 key_function=get_ccds_accession)

    # parse gnomad transcript-level variants
    print('Reading gnomAD variant database ...')
    with open(configs['gnomad_variants'], 'rt') as variant_handle:
        # transcript_variants will be a dict of dicts where major version
        # ENSEMBL transcript IDs are the first level keys and "ccds", "ensp",
        # "swissprot", "variants" are the second level keys. The value of each
        # second-level key is a Python list.
        transcript_variants = json.load(variant_handle)


    output_dir = os.path.abspath(configs['output_dir'])


    # compute the MRT scores for each transcript
    rate_vs_count = [] 
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
            if any([x not in {'A', 'T', 'C', 'G'} for x in set(transcript_cds)]):
                print('Invalid CDS! Skipped.')
                print(transcript_cds)
                continue

            # calculate expected counts for each codon
            try:
                codon_mutation_rates = get_codon_mutation_rates(transcript_cds)
            except ValueError:
                continue

            # determine the total rate for the current transcript
            synonymous_rate = 0
            nonsynonymous_rate = 0
            for rates in codon_mutation_rates:
                synonymous_rate += rates[0]
                nonsynonymous_rate += rates[1]

            # get all variants of this transcript reported in gnomAD
            try:
                variants = transcript_variants[transcript]['variants']
            except KeyError:
                logging.critical('No variants found in %s in gnomAD', transcript)
                logging.critical('%s was skipped ...', transcript)
                continue

            # tabulate variants at each site
            synonymous_count, nonsynonymous_count = count_variants(variants)

            print(transcript, len(transcript_cds) // 3, 
                '{:.3e}'.format(synonymous_rate), 
                synonymous_count,
                '{:.3e}'.format(nonsynonymous_rate), 
                nonsynonymous_count)
            rate_vs_count.append([
                transcript, 
                len(transcript_cds) // 3, 
                '{:.3e}'.format(synonymous_rate), 
                synonymous_count,
                '{:.3e}'.format(nonsynonymous_rate), 
                nonsynonymous_count
            ]) 

        with open(file=os.path.join(output_dir, args.output), mode='wt') as opf:
            header = ['TRANSCRIPT', 'LENGTH', 'SYNONYMOUS_RATE', 'SYNONYMOUS_COUNT', 
                'MISSENSE_RATE', 'MISSENSE_COUNT']
            csv_writer = csv.writer(opf, delimiter='\t')
            csv_writer.writerow(header)
            csv_writer.writerows(rate_vs_count)


if __name__ == '__main__':
    main()
