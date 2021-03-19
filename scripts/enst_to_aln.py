#!/usr/bin/env python3

from Bio import SeqIO
from Bio import BiopythonWarning
import gzip
import warnings
import logging
from argparse import ArgumentParser
from collections import defaultdict

# ignore warnings from Biopython
warnings.simplefilter('ignore', BiopythonWarning)


def parse_cmd_args():
    """

    Returns
    -------

    """
    parser = ArgumentParser(description='''To be added ...''')
    parser.add_argument(
        '-i', '--input', dest='transcripts', type=str, required=True,
        help='''A lists of Ensembl transcript IDs for which the alignment is to
        be extracted.'''
    )
    parser.add_argument(
        '-s', '--seq-file', dest='seq_file', type=str, required=True,
        help='''UCSC 100/30-way sequence alignment file.'''
    )
    parser.add_argument(
        '-w', '--way', dest='way', type=int, required=True, choices=[30, 100],
        help='''100- or 30-way alignment. This parameter is needed because
        USCS splits the 100-way alignments into exons. The script needs to
        knit the exon sequences together for each transcript.'''
    )
    parser.add_argument(
        '--seq-type', dest='seq_type', type=str, default='nuc',
        choices=['nuc', 'aa'], help='''The type of input sequence.'''
    )
    parser.add_argument(
        '--log', dest='log', type=str, required=False, default='enst_to_aln.log',
        help='''Log file to store messages from running this script.'''
    )
    return parser.parse_args()


def main():
    """

    Returns
    -------

    """
    # parse command-line arguments
    args = parse_cmd_args()

    # configure the logger
    logging.basicConfig(
        filename=args.log,
        level=logging.CRITICAL,
        filemode='w',
        format='%(levelname)s:%(asctime)s:%(message)s'
    )

    # read sequence file
    alignments = {}

    # read in UCSC multiple sequence alignments
    print('Parsing sequence alignments ...')
    with gzip.open(args.seq_file, 'rt') as ip_handle:
        species_seqs = defaultdict(str)
        for seq_record in SeqIO.parse(ip_handle, format='fasta'):
            parts = seq_record.description.split()
            id_fields = parts[0].split('_')
            enst_id_full, species = id_fields[:2]
            enst_id_major = enst_id_full.split('.')[0]
            if args.way == 30:
                species_seqs[species] = seq_record.seq

                # hit the sequence of the last species
                if species == 'dasNov3':
                    alignments[enst_id_major] = species_seqs
                    # reset for the next transcript
                    species_seqs = defaultdict(str)

            else:
                cur_exon, last_exon = [int(x) for x in id_fields[2:4]]
                if cur_exon <= last_exon:
                    species_seqs[species] += seq_record.seq

                if cur_exon == last_exon and species == 'petMar2':
                    alignments[enst_id_major] = species_seqs
                    # reset for the next transcript
                    species_seqs = defaultdict(str)

    # print the number of parsed alignments
    print('Parsed {} alignments from {}.'.format(len(alignments), args.seq_file))

    # read in all the Ensembl transcript IDs
    with open(args.transcripts, 'rt') as ip_handle:
        input_transcript_ids = [line.strip() for line in ip_handle]
    print('Read {} transcripts in.'.format(len(input_transcript_ids)))

    # write the multiple sequence alignment for each transcript in FASTA format
    for enst_id in input_transcript_ids:
        try:
            enst_alignment = alignments[enst_id]
        except KeyError:
            logging.critical(
                'No alignment found for {}.'.format(enst_id)
            )
            continue

        aln_file = enst_id + '_aln.fasta'
        with open(aln_file, 'wt') as op_handle:
            if args.seq_type == 'nuc':
                offset = 3
            else:
                offset = 1
            for species, seq in enst_alignment.items():
                op_handle.write('>{}\n{}\n'.format(species, seq[:-offset]))


if __name__ == '__main__':
    main()
