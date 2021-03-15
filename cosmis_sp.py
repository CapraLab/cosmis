#!/usr/bin/env python3

"""
    @summary: This is the single-protein version of the COSMIS method, hence its
    name 'cosmis_sp.py'. It was designed to work with proteins for which the
    COSMIS scores were not already computed, because the protein structure
    database used wasn't up to date. However, if there exists a structure or
    homology model for the protein of interest, its COSMIS scores can still be
    computed using this script.
    This script assumes that residues in the PDB file are numbered according to
    the amino acid sequence in UniProt.
    @author: Bian Li
    @contact: bian.li@vanderbilt.edu
    @change: Last modified 2/12/2021.

"""

import sys, csv
import gzip
import json
from collections import defaultdict
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.SeqUtils import seq1
from Bio.PDB import PDBParser, is_aa
from cosmis.utils import pdb_utils, seq_utils


def parse_cmd():
    """
    Specifies command-line flags and parses command-line options.

    Returns
    -------
    ArgumentParser
        An object of type ArgumentParser containing the parsed command-line
        arguments.

    """
    parser = ArgumentParser()
    parser.add_argument('-c', '--config', dest='config', required=True,
                        type=str, help='A JSON file specifying options.')
    parser.add_argument('-t', '--transcript', dest='transcript', required=True,
                        type=str, help='''Ensembl transcript ID of the 
                        transcript for which to compute a COSMIS profile.''')
    parser.add_argument('-p', '--pdb', dest='pdb_file', required=True,
                        type=str, help='''PDB file containing the structure
                        of the protein structure of the given transcript.''')
    parser.add_argument('-o', '--output', dest='output_file', required=True,
                        type=str, help='''Output file to store the COSMIS scores
                        of the protein.''')
    parser.add_argument('--multimer', dest='multimer', action='store_true',
                        default=False, help='Is the input PDB file a multimer?')
    parser.add_argument('--chain', dest='pdb_chain', default='A', type=str,
                        help='Chain ID of the subunit in the PBD file.')
    parser.add_argument('-w', '--overwrite', dest='overwrite', required=False,
                        action='store_true', help='''Whether to overwrite 
                        already computed MTR3D scores.''')
    parser.add_argument('-v', '--verbose', dest='verbose', required=False,
                        action='store_true', help='''Whether to output verbose
                        data: number of contacting residues and number of 
                        missense and synonymous variants in the neighborhood
                        of the mutation site.''')
    return parser.parse_args()


def get_ensembl_accession(record):
    """
    Convert a FASTA file record into an accession ID.

    Parameters
    ----------
    record : str
        Record ID is of the format: ">CCDS2.2|Hs109|chr1"

    Returns
    -------
    str
        Accession ID that can be used as dictionary key.

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


def get_transcript_pep_seq(enst_id, ensp_id, pep_dict):
    """

    Parameters
    ----------
    enst_id : str

    ensp_id : str

    pep_dict : dict

    Returns
    -------

    """
    try:
        transcript_pep = pep_dict[ensp_id].seq
    except KeyError:
        print('%s not found in given database' % ensp_id)
        print('%s was skipped ...' % enst_id)
        return None
    return transcript_pep


def parse_config(config):
    """
    Parses the configuration file for computing COSMIS scores into a Python
    dictionary.

    Parameters
    ----------
    config : json
        Configuration file in JSON format.

    Returns
    -------
    dict
        A Python dictionary.

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
    dict
        A dictionary where the key is amino acid position and the value is
        the number of variants at this position. One dictionary for missense
        variants and one dictionary for synonymous variants.

    """
    #
    missense_counts = defaultdict(int)
    synonymous_counts = defaultdict(int)
    for variant in variants:
        vv, ac, an = variant
        w = vv[0]  # wild-type amino acid
        v = vv[-1]  # mutant amino acid
        pos = vv[1:-1]  # position in the protein sequence
        # only consider rare variants
        if int(ac) / int(an) > 0.001:
            continue
        if w != v:  # missense variant
            missense_counts[int(pos)] += 1
        else:  # synonymous variant
            synonymous_counts[int(pos)] += 1
    return missense_counts, synonymous_counts


def parse_ensembl_cds(file=None):
    """
    
    Parameters
    ----------
    file

    Returns
    -------

    """
    print('Parsing ENSEMBL concensus coding sequence database ...')
    with gzip.open(file, mode='rt') as cds_handle:
        ensembl_cds_dict = SeqIO.to_dict(
            SeqIO.parse(cds_handle, format='fasta'),
            key_function=get_ensembl_accession
        )
    return ensembl_cds_dict


def parse_ccds(file=None):
    """

    Parameters
    ----------
    file

    Returns
    -------

    """
    print('Parsing NCBI CCDS database ...')
    with gzip.open(file, mode='rt') as ccds_handle:
        ccds_dict = SeqIO.to_dict(
            SeqIO.parse(ccds_handle, format='fasta'),
            key_function=get_ccds_accession
        )
    return ccds_dict


def parse_ensembl_pep(file=None):
    """

    Parameters
    ----------
    file

    Returns
    -------

    """
    print('Parsing Ensembl protein sequence database ...')
    with gzip.open(file, mode='rt') as pep_handle:
        pep_dict = SeqIO.to_dict(
            SeqIO.parse(pep_handle, format='fasta'),
            key_function=get_ensembl_accession
        )
    return pep_dict


def main():
    # parse command-line arguments
    args = parse_cmd()

    # parse configuration file
    configs = parse_config(args.config)

    # ENSEMBL cds
    ensembl_cds_dict = parse_ensembl_cds(configs['ensembl_cds'])
    ccds_dict = parse_ccds(configs['ccds_cds'])
    pep_dict = parse_ensembl_pep(configs['ensembl_pep'])

    # parse gnomad transcript-level variants
    print('Reading gnomAD variant database ...')
    with open(configs['gnomad_variants'], 'rt') as variant_handle:
        # transcript_variants will be a dict of dicts where major version
        # ENSEMBL transcript IDs are the first level keys and "ccds", "ensp",
        # "swissprot", "variants" are the second level keys. The value of each
        # second-level key is a Python list.
        transcript_variants = json.load(variant_handle)

    # compute the MRT scores
    pdb_file = args.pdb_file
    pdb_chain = args.pdb_chain
    transcript = args.transcript

    # get the amino acid sequence of the transcript
    try:
        # Ensembl peptide ID for the transcript
        ensp_id = transcript_variants[transcript]['ensp'][0]
    except KeyError:
        print(
            'Transcript %s not found in %s', 
            transcript, 
            configs['gnomad_variants']
        )
        sys.exit(1)

    # get the peptide sequence from peptide sequence database
    cosmis_scores = []
    transcript_pep_seq = get_transcript_pep_seq(
        transcript, ensp_id, pep_dict
    )

    # get all variants of this transcript reported in gnomAD
    try:
        variants = transcript_variants[transcript]['variants']
    except KeyError:
        print('No variants found for %s in gnomAD', transcript)
        sys.exit(1)

    # get the coding sequence of the transcript
    try:
        transcript_cds = ensembl_cds_dict[transcript]
    except KeyError:
        print(
            '''No CDS found in Ensembl CDS database! 
            Looking for it in the CCDS database ...'''
        )
        transcript_cds = None

    if transcript_cds is None:
        try:
            ccds_id = transcript_variants[transcript]['ccds'][0]
            transcript_cds = ccds_dict[ccds_id]
        except KeyError:
            print('ERROR: No CDS found in CCDS database!')
            sys.exit(1)

    # check that the CDS does not contain invalid nucleotides
    if seq_utils.is_valid_cds(transcript_cds):
        print('ERROR: invalid CDS for', transcript_cds)
        sys.exit(1)

    uniprot_ids = transcript_variants[transcript]['swissprot']
    if len(uniprot_ids) > 1:
        print(
            'ERROR: more than one UniProt IDs found for transcript '
            '%s: %s', transcript, ','.join(uniprot_ids)
        )
        sys.exit(1)

    # get the UniProt ID
    uniprot_id = uniprot_ids[0]

    # print message
    print(
        'Estimating COSMIS scores for:',
        transcript, 
        ensp_id,
        pdb_file
    )

    pdb_parser = PDBParser(PERMISSIVE=1)
    structure = pdb_parser.get_structure(id='NA', file=pdb_file)
    chain = structure[0][pdb_chain]
    if args.multimer:
        all_aa_residues = [aa for aa in structure[0].get_residues() if is_aa(aa)]
    else:
        all_aa_residues = [aa for aa in chain.get_residues() if is_aa(aa)]
    all_contacts = pdb_utils.search_for_all_contacts(all_aa_residues, radius=8)

    # get mutation rates and expected counts for each codon
    codon_mutation_rates = seq_utils.get_codon_mutation_rates(transcript_cds.seq)
    all_cds_ns_counts = seq_utils.count_poss_ns_variants(transcript_cds.seq)
    
    # tabulate variants at each site
    # missense_counts and synonymous_counts are dictionary that maps
    # amino acid positions to variant counts
    missense_counts, synonymous_counts = count_variants(variants)

    # index all contacts by residue ID
    indexed_contacts = defaultdict(list)
    for c in all_contacts:
        indexed_contacts[c.get_res_a()].append(c.get_res_b())
        indexed_contacts[c.get_res_b()].append(c.get_res_a())

    for seq_pos, seq_aa in enumerate(transcript_pep_seq, start=1):
        # check that the amino acid in ENSP sequence matches 
        # that in the PDB structure
        try:
            res = chain[seq_pos]
        except KeyError:
            print(
                'Residue %s not found in chain %s in PDB file: %s' % 
                (seq_pos, pdb_chain, pdb_file)
            )
            continue
        pdb_aa = seq1(res.get_resname())
        if seq_aa != pdb_aa:
            print('Residue in ENSP did not match that in PDB at', seq_pos)
            continue

        contact_res = indexed_contacts[res]
        contacts_pdb_pos = [r.get_id()[1] for r in contact_res]
        
        seq_seps = ';'.join(
            str(x) for x in [i - seq_pos for i in contacts_pdb_pos]
        )

        total_missense_obs = missense_counts.setdefault(seq_pos, 0)
        total_synonymous_obs = synonymous_counts.setdefault(seq_pos, 0)
        total_missense_poss = all_cds_ns_counts[seq_pos - 1][0]
        total_synonyms_poss = all_cds_ns_counts[seq_pos - 1][1]
        total_synonymous_rate = codon_mutation_rates[seq_pos - 1][0]
        total_missense_rate = codon_mutation_rates[seq_pos - 1][1]
        for j in contacts_pdb_pos:
            # count the total # observed variants in contacting residues
            total_missense_obs += missense_counts.setdefault(j, 0)
            total_synonymous_obs += synonymous_counts.setdefault(j, 0)

            # count the total # expected variants
            total_missense_poss += all_cds_ns_counts[j - 1][0]
            total_synonyms_poss += all_cds_ns_counts[j - 1][1]
            total_synonymous_rate += codon_mutation_rates[j - 1][0]
            total_missense_rate += codon_mutation_rates[j - 1][1]

        # compute the fraction of expected missense variants
        cosmis_scores.append(
            [
                transcript, ensp_id, uniprot_id, seq_pos, seq_aa, seq_seps,
                len(contacts_pdb_pos), total_synonyms_poss, total_missense_poss,
                '%.3e' % total_synonymous_rate, total_synonymous_obs,
                '%.3e' % total_missense_rate, total_missense_obs
            ]
        )

    with open(file=args.output_file, mode='wt') as opf:
        header = [
            'enst_id', 'ensp_id', 'uniprot_id', 'ensp_pos', 'ensp_aa',
            'seq_separations', 'num_contacts', 'synonymous_poss',
            'missense_poss', 'synonymous_rate', 'synonymous_obs',
            'missense_rate', 'missense_obs'
        ]
        csv_writer = csv.writer(opf, delimiter='\t')
        csv_writer.writerow(header)
        csv_writer.writerows(cosmis_scores)


if __name__ == '__main__':
    main()
