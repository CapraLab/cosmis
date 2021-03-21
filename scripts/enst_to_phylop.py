#!/usr/bin/env python3

import json
import gzip
from Bio import SeqIO
from argparse import ArgumentParser
from collections import defaultdict
from cosmis.utils import seq_utils


def parse_cmd_args():
    """
    Parse command-line arguments.

    Returns
    -------
    Parsed command-line arguments.

    """
    parser = ArgumentParser(description='')
    parser.add_argument(
        '-c', '--coordinates', dest='coordinates', type=str, required=True,
        help='''Input file containing the genomic coordinates of each of the
        GENCODE annotated transcript coding sequences. One per line, 
        preferably in BED format.'''
    )
    parser.add_argument(
        '-f', '--format', dest='format', type=str, required=False, default='bed',
        choices=['bed', 'gencode'],
        help='''Format of the input coordinate file. This is used to 
        determine whether the coordinates are 0-based or 1-based.'''
    )
    parser.add_argument(
        '-p', '--phylop', dest='phylop', type=str, required=True,
        help='Basewise phyloP scores indexed by 1-based genomic coordinates.'
    )
    parser.add_argument(
        '-s', '--sequence-file', dest='sequence_file', type=str, required=True,
        help='''A gzip file containing the coding nucleotide sequences of all 
        GENCODE annotated transcripts.'''
    )
    parser.add_argument(
        '-v', '--verbose', dest='verbose', action='store_true',
        help='''Whether to output detailed messages while the script is running.'''
    )
    parser.add_argument(
        '-o', '--output', dest='output', type=str, required=False,
        default='out.json', help='''Output file.'''
    )
    # do any necessary checks on command-line arguments here
    return parser.parse_args()


def main():
    # parse all command-line arguments
    cmd_args = parse_cmd_args()

    coords_dict = defaultdict(list)
    chrom_dict = {}
    # parse coordinates into a dictionary
    with open(cmd_args.coordinates, 'rt') as ipf:
        for line in ipf:
            chrom, start, end, strand, enst = line.strip().split()
            start = int(start)
            end = int(end)
            enst_id = enst.split('.')[0]
            # add the genomic coordinates for this CDS
            if cmd_args.format == 'bed':
                start += 1
            # using 1-start, fully-closed coordinate system
            if chrom_dict.get(enst_id) is None:
                chrom_dict[enst_id] = chrom
            if strand == '-':
                start, end = end, start
                coords_dict[enst_id].extend(range(start, end - 1, -1))
            else:
                coords_dict[enst_id].extend(range(start, end + 1))

    print('Extracted coordinates for {} transcripts.'.format(len(coords_dict)))

    # needs to use the gzip module to help read compressed sequences
    seqs_dict = {}
    with gzip.open(cmd_args.sequence_file, 'rt') as handle:
        for seq_record in SeqIO.parse(handle, format='fasta'):
            # find the CDS field
            id_fields = seq_record.id.split('|')
            for f in id_fields:
                if f.startswith('CDS:'):
                    cds_field = f
                    break
            cds_start, cds_end = [
                int(x) for x in cds_field.split(':')[1].split('-')
            ]
            cds = seq_record[(cds_start-1):cds_end]
            # translate CDS to amino acid sequence
            if not seq_utils.is_valid_cds(cds.seq):
                continue
            aa_seq = seq_utils.translate(cds)
            enst_id = id_fields[0].split('.')[0]
            seqs_dict[enst_id] = (cds, aa_seq)

    print('Extracted {} canonical CDS records.'.format(len(seqs_dict)))

    phylop_dict = defaultdict(dict)
    with gzip.open(cmd_args.phylop, 'rt') as ipf:
        for line in ipf:
            chrom, start, end, score = line.strip().split()
            if not chrom.startswith('chr'):
                chrom = 'chr' + chrom
            if int(end) - int(start) == 1:
                phylop_dict[chrom][int(end)] = float(score)
            # this line represents multiple bases collapsed into an interval
            else:
                for i in range(int(start) + 1, int(end) + 1):
                    phylop_dict[chrom][i] = float(score)

    print('Extracted {} phyloP scores.'.format(len(phylop_dict)))

    # map CDS to genomic coordinates
    all_cds_to_phylop = defaultdict(dict)
    for enst_id, (cds, aa_seq) in seqs_dict.items():
        # get genomic coordinates
        try:
            coords = coords_dict[enst_id]
            chrom = chrom_dict[enst_id]
        except KeyError:
            if cmd_args.verbose:
                print('No coordinates found for {}.'.format(enst_id))
                print('Please double check your input files.')
            continue

        # add chrom info
        all_cds_to_phylop[enst_id]['chrom'] = chrom

        if len(cds) - 3 != len(coords):
            if cmd_args.verbose:
                print('{} CDS length does not match the number of coordinates.'.format(enst_id))
                print('CDS length: {}, # coordinates: {}'.format(len(cds), len(coords)))
            continue

        cds_to_phylop = []
        for i, a in enumerate(aa_seq):
            codon = str(cds[i*3:(i*3+3)].seq)
            codon_coords = coords[i*3:(i*3+3)]
            try:
                codon_phylop_scores = [
                    phylop_dict.get(chrom).get(int(x), '') for x in codon_coords
                ]
            except AttributeError:
                print('No record found for {}:{}.'.format(chrom, codon_coords))
                codon_phylop_scores = ['', '', '']
            cds_to_phylop.append(
                tuple([a, codon, codon_coords, codon_phylop_scores])
            )
        all_cds_to_phylop[enst_id]['gerp'] = cds_to_phylop

    print('Now writing mapping to {}'.format(cmd_args.output))
    with open(cmd_args.output, 'wt') as opf:
        json.dump(all_cds_to_phylop, fp=opf, indent=4)


if __name__ == '__main__':
    main()
