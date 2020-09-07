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
        '-g', '--gene-list', dest='gene_list', required=True, type=str,
         help='''A list of Ensembl gene IDs for which the canonical 
         transcripts are to be retrieved.'''       
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
    
    id_fields = seq_record.id.split('|')
    
    # find the CDS field
    for f in id_fields:
        if f.startswith('CDS'):
            cds_field = f
            break
    try:
        cds_start, cds_end = [
            int(x) for x in cds_field.split(':')[1].split('-')
        ]
    except (ValueError, IndexError):
        print('ERROR in retrieving CDS starting and ending positions.')
        print(cds_field)
        return False
                
    # the following conditionals skip noncanonical transcripts 
    if not seq_record.seq[cds_start-1:cds_start+2] == 'ATG':
        print('The CDS of', id_fields[0], 'does not start with ATG.')
        return False
    if seq_record.seq[cds_end-3:cds_end] not in stop_codons:
        print('The CDS of', id_fields[0], 'does not have a stop codon.')
        return False
    if not (cds_end - cds_start + 1) % 3 == 0:
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
            # keep cds of valid transcripts
            if is_valid_transcript(seq_record):
                # find the CDS field
                id_fields = seq_record.id.split('|')
                for f in id_fields:
                    if f.startswith('CDS'):
                        cds_field = f
                        break
                cds_start, cds_end = [
                    int(x) for x in cds_field.split(':')[1].split('-')
                ]
                cds = seq_record[cds_start-1:cds_end]
                ensembl_gene_id = id_fields[1].split('.')[0]
                gene_to_transcripts[ensembl_gene_id].append(cds)
            
    # read in the list of gene IDs
    with open(args.gene_list, 'rt') as ipf:
        gene_list = [l.strip() for l in ipf]
        
    # now retrieve canonical transcript for each gene
    canonical_transcripts = []
    
    for gene in gene_list:
        # skip genes for which no transcripts are available
        transcripts = gene_to_transcripts[gene]
        if not transcripts:
            print('No transcript found for', gene)
            continue
        
        # canonical is defined as the longest transcript
        canonical_transcript = max(transcripts, key=lambda t: len(t))
        canonical_transcripts.append(canonical_transcript)
        
    # write canonical transcripts to disk file
    with gzip.open(args.output, 'wt') as handle:
        SeqIO.write(canonical_transcripts, handle, 'fasta')
        

if __name__ == '__main__':
    main()
