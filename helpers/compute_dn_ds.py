#!/usr/bin/env python3


from Bio import SeqIO
import gzip
from argparse import ArgumentParser
from collections import defaultdict
import os
from Bio.Phylo.PAML import codeml


def parse_cmd_args():
    """

    Returns
    -------

    """
    parser = ArgumentParser(description='''To be added ...''')
    parser.add_argument(
        '-i', '--input', dest='transcripts', type=str, required=True,
        help='''A lists of Ensembl transcript IDs for each the dN/dS is to be
        computed.'''
    )
    parser.add_argument(
        '-s', '--seq-file', dest='seq_file', type=str, required=True,
        help='''UCSC 100-way sequence alignment file.'''
    )
    parser.add_argument(
        '-t', '--tree', dest='tree', type=str, required=True,
        help='''UCSC 100-way phylogenetic tree file.'''
    )
    parser.add_argument(
        '-c', '--control-file', dest='control_file', type=str,
        default='codeml.ctl', help='''CODEML control file.'''
    )
    parser.add_argument(
        '-b', '--codeml', dest='codeml', type=str, default='codeml',
        help='''Path to the binary of CODEML.'''
    )
    parser.add_argument(
        '--seq-type', dest='seq_type', type=str, default='nuc',
        choices=['nuc', 'aa'], help='''The type of input sequence.'''
    )
    return parser.parse_args()


def main():
    """

    Returns
    -------

    """
    # parse command-line arguments
    args = parse_cmd_args()

    # read sequence file
    alignments = defaultdict(dict)

    #
    print('Parsing sequence alignments ...')
    with gzip.open(args.seq_file, 'rt') as ip_handle:
        species_seqs = defaultdict(str)
        for seq_record in SeqIO.parse(ip_handle, format='fasta'):
            parts = seq_record.description.split()
            enst_id_full, species = parts[0].split('_')
            enst_id_major = enst_id_full.split('.')[0]
            species_seqs[species] = seq_record.seq

            # hit the last exon of the last species
            if species == 'dasNov3':
                # last exon for this transcript
                alignments[enst_id_major]['length'] = int(parts[1])
                alignments[enst_id_major]['alignment'] = species_seqs
                species_seqs = defaultdict(str)

    #
    with open(args.transcripts, 'rt') as ip_handle:
        input_transcript_ids = [l.strip() for l in ip_handle]
    #
    print('Read {} transcripts in.'.format(len(input_transcript_ids)))

    # run codeml for each transcript
    for transcript in input_transcript_ids:
        with open('seqfile.txt', 'wt') as op_handle:
            if args.seq_type == 'nuc':
                offset = 3
            else:
                offset = 1
            op_handle.write(' 30 {}\n'.format(
                alignments[transcript]['length'] - offset
            ))
            for species, seq in alignments[transcript]['alignment'].items():
                op_handle.write('{}\n{}\n'.format(species, seq[:-offset]))

        cml = codeml.Codeml(
            alignment='seqfile.txt',
            tree=args.tree,
            out_file='results.txt'
        )
        cml.read_ctl_file(args.control_file)

        #
        print('Now running codeml for transcript: {}'.format(transcript))
        print('========================================================')
        cml.run(verbose=True, command=args.codeml)

        # rename files
        if os.path.exists('rst'):
            os.rename(src='rst', dst=transcript + '_rst')
        else:
            print('"rst" file does not exist.')
        if os.path.exists('rst1'):
            os.rename(src='rst1', dst=transcript + '_rst1')
        else:
            print('"rst1" file does not exist.')
        if os.path.exists('results.txt'):
            os.rename(src='results.txt', dst=transcript + '_results.txt')
        else:
            print('"results.txt" file does not exist.')
        if os.path.exists('seqfile.txt'):
            os.rename(src='seqfile.txt', dst=transcript + '_seqfile.txt')
        else:
            print('"seqfile.txt" file does not exist.')


if __name__ == '__main__':
    main()
