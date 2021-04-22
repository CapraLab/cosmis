#!/usr/bin/env python3

"""
For a given list of Ensembl gene IDs, this script extract the canonical
transcript for each of the genes.

@author: Bian Li
@email: bian.li@vanderbilt.edu
"""

from Bio import SeqIO
import gzip
from argparse import ArgumentParser
from collections import defaultdict


def parse_cmd_args():
    """
    Specifies arguments that users can/should give on the command line.

    Returns
    -------
    An object of type ArgumentParser.

    """
    parser = ArgumentParser(
        description='''For a given list of Ensembl gene IDs, this script
        extract the canonical transcript for each of the genes.'''
    )
    parser.add_argument(
        '-i', '--input', dest='input', required=True, type=str,
        help='''A gzip compressed file containing GENCODE annotated protein-
        coding transcript sequences.'''
    )
    parser.add_argument(
        '-o', '--output', dest='output', required=True, type=str,
        help='Name of the disk file to store the processed sequences.'
    )
    parser.add_argument(
        '-f', '--format', dest='format', required=False, type=str,
        default='fasta', help='Format of the sequence file.'
    )

    # return all command-line arguments
    return parser.parse_args()


def is_valid_transcript(seq_record):
    """
    Determines whether a given sequence is a valid transcript sequence or not.
    A valid transcript sequence is defined as one that starts with ATG, ends
    with a stop codon, that whose length is a multiple of three.

    Parameters
    ----------
    seq_record : SeqRecord
        An object of type SeqRecord.

    Returns
    -------
    bool
        True if the given sequence is a valid transcript sequence else False.

    """
    # stop codons
    stop_codons = {'TAG', 'TAA', 'TGA'}

    # get all ID fields
    id_fields = seq_record.description.split(' ')

    # the following conditionals skip noncanonical transcripts
    if seq_record.seq[:3] != 'ATG':
        print('The CDS of', id_fields[0], 'does not start with ATG.')
        return False
    if seq_record.seq[-3:] not in stop_codons:
        print('The CDS of', id_fields[0], 'does not have a stop codon.')
        return False
    if len(seq_record) % 3 != 0:
        print('The CDS of', id_fields[0], 'is not a multiple of 3.')
        return False

    return True


def main():
    # return all command-line arguments
    args = parse_cmd_args()

    #
    gene_to_transcripts = defaultdict(list)

    # needs to use the gzip module to help read compressed sequences
    with gzip.open(args.input, 'rt') as handle:
        for seq_record in SeqIO.parse(handle, args.format):
            id_fields = seq_record.description.split(' ')
            # keep cds of valid transcripts
            if is_valid_transcript(seq_record):
                try:
                    ensembl_gene_id = id_fields[3].split(':')[1].split('.')[0]
                except IndexError:
                    print('BAD fasta header!', seq_record.description)
                    continue
                gene_to_transcripts[ensembl_gene_id].append(seq_record)

    # now retrieve canonical transcript for each gene
    canonical_transcripts = []

    for gene in gene_to_transcripts.keys():
        # canonical is defined as the longest transcript
        transcripts = gene_to_transcripts[gene]
        canonical_transcript = max(transcripts, key=lambda t: len(t))
        canonical_transcripts.append(canonical_transcript)

    # write canonical transcripts to disk file
    with gzip.open(args.output, 'wt') as handle:
        SeqIO.write(canonical_transcripts, handle, 'fasta')


if __name__ == '__main__':
    main()
