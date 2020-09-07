#!/usr/bin/env python3

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
from Bio import SeqIO
from Bio.SeqUtils import seq1
from Bio.PDB import PDBParser, NeighborSearch, is_aa, PDBList

from Bio import BiopythonWarning
import warnings
warnings.simplefilter('ignore', BiopythonWarning)

GENETIC_CODE = {
    'TTT': 'F',
    'CTT': 'L',
    'ATT': 'I',
    'GTT': 'V',
    'TTC': 'F',
    'CTC': 'L',
    'ATC': 'I',
    'GTC': 'V',
    'TTA': 'L',
    'CTA': 'L',
    'ATA': 'I',
    'GTA': 'V',
    'TTG': 'L',
    'CTG': 'L',
    'ATG': 'M',
    'GTG': 'V',
    'TCT': 'S',
    'CCT': 'P',
    'ACT': 'T',
    'GCT': 'A',
    'TCC': 'S',
    'CCC': 'P',
    'ACC': 'T',
    'GCC': 'A',
    'TCA': 'S',
    'CCA': 'P',
    'ACA': 'T',
    'GCA': 'A',
    'TCG': 'S',
    'CCG': 'P',
    'ACG': 'T',
    'GCG': 'A',
    'TAT': 'Y',
    'CAT': 'H',
    'AAT': 'N',
    'GAT': 'D',
    'TAC': 'Y',
    'CAC': 'H',
    'AAC': 'N',
    'GAC': 'D',
    'TAA': 'Stop',
    'CAA': 'Q',
    'AAA': 'K',
    'GAA': 'E',
    'TAG': 'Stop',
    'CAG': 'Q',
    'AAG': 'K',
    'GAG': 'E',
    'TGT': 'C',
    'CGT': 'R',
    'AGT': 'S',
    'GGT': 'G',
    'TGC': 'C',
    'CGC': 'R',
    'AGC': 'S',
    'GGC': 'G',
    'TGA': 'Stop',
    'CGA': 'R',
    'AGA': 'R',
    'GGA': 'G',
    'TGG': 'W',
    'CGG': 'R',
    'AGG': 'R',
    'GGG': 'G'
}


def count_codon_ns(codon):
    """

    Parameters
    ----------
    codon : str
        A three-letter genetic codon.

    Returns
    -------
    tuple
        A pair of # missense variants vs # synonymous variants.

    """
    nucleotides = {'A', 'T', 'C', 'G'}
    aa = GENETIC_CODE[codon]
    missense = 0
    synonymous = 0
    for i, x in enumerate(codon):
        for n in nucleotides - {x}:
            new_codon = codon[:i] + n + codon[i+1:]
            new_aa = GENETIC_CODE[new_codon]
            if new_aa == aa:
                synonymous += 1
            elif new_aa != aa and new_aa != 'Stop':
                missense += 1
    return missense, synonymous


def count_cds_ns(cds):
    """

    Parameters
    ----------
    cds : str
        Coding sequence.

    Returns
    -------
    list
        A pair of # missense variants vs # synonymous variants for each
        codon in the coding sequence.

    """
    cds_ns = []
    for i in range(0, len(cds), 3):
        cds_ns.append(count_codon_ns(cds[i:i+3]))
    return cds_ns


def search_for_all_contacts(residues, radius=6):
    """
    Search for all contacts in the given set of residues based on
    distances between heavy atoms (atoms other than hydrogen atoms).

    Parameters
    ----------
    residues
    radius

    Returns
    -------

    """
    atom_list = []
    for r in residues:
        if r.get_resname() == 'GLY':
            try:
                atom_list.append(r['CA'])
            except KeyError:
                print('No CA atom found for GLY:', r, 'skipped ...')
                continue
        else:
            atom_list += [a for a in r.get_atoms() if a.get_name() 
                          not in BACKBONE_ATOMS]
        # atom_list += [a for a in r.get_atoms()]
    ns = NeighborSearch(atom_list)
    all_contacts = [Contact(res_a=c[0], res_b=c[1]) for c in ns.search_all(radius, level='R')]
    return all_contacts


def compute_mtr1d(pos, ns_counts, syn_counts, expected_counts, window=31):
    """

    Parameters
    ----------
    pos : int
        Sequence position.
    ns_counts : dict

    syn_counts : dict

    expected_counts : list

    window : int
        Window size.

    Returns
    -------
    float


    """
    total_ns_obs = 0
    total_syn_obs = 0
    total_ns_exp = 0
    total_syn_exp = 0
    # handle the special case of the first 15 residues
    if pos < 16:
        for i in range(31):
            try:
                total_ns_obs += ns_counts[i]
            except KeyError:
                total_ns_obs += 0
            try:
                total_syn_obs += syn_counts[i]
            except KeyError:
                total_syn_obs += 0
            total_ns_exp += expected_counts[i][0]
            total_syn_exp += expected_counts[i][1]
    # handle the special case of the last 15 residues
    elif pos > len(expected_counts) - 15:
        for i in range(len(expected_counts) - 31, len(expected_counts)):
            try:
                total_ns_obs += ns_counts[i]
            except KeyError:
                total_ns_obs += 0
            try:
                total_syn_obs += syn_counts[i]
            except KeyError:
                total_syn_obs += 0
            total_ns_exp += expected_counts[i][0]
            total_syn_exp += expected_counts[i][1]
    else:
        for i in range(pos - window // 2, pos + window // 2 + 1):
            try:
                total_ns_obs += ns_counts[i]
            except KeyError:
                total_ns_obs += 0
            try:
                total_syn_obs += syn_counts[i]
            except KeyError:
                total_syn_obs += 0
            if i > 1 or i < len(expected_counts):
                total_ns_exp += expected_counts[i - 1][0]
                total_syn_exp += expected_counts[i - 1][1]
    try:
        mtr1d = (total_ns_obs / (total_ns_obs + total_syn_obs)) / (total_ns_exp / (total_ns_exp + total_syn_exp))
    except ZeroDivisionError:
        mtr1d = 0

    return total_ns_obs, total_syn_obs, mtr1d


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
    parser.add_argument('-p', '--pdb', dest='pdb', type=str, required=False,
                        help='PDB file representing the 3D structure of the '
                             'protein')
    parser.add_argument('-w', '--overwrite', dest='overwrite', required=False,
                        action='store_true', help='Whether to overwrite already '
                        'computed MTR3D scores.')
    parser.add_argument('--mtr1d', dest='mtr1d', required=False, default=False,
                        action='store_true', help='If specified, computeds MTR '
                        'score implemented according to Traynelis et al., Genome '
                        'Research, 2017')
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
    # if len(ensp_id) > 1:
    #     logging.critical(
    #         'More than one ENSP IDs found for transcript %s: %s',
    #         enst_id, ensp_id
    #     )
    #     sys.exit(1)
    try:
        transcript_pep = pep_dict[ensp_id].seq
    except KeyError:
        logging.critical('%s not found in given database', ensp_id)
        logging.critical('%s was skipped ...', enst_id)
        print('%s not found in given database' % ensp_id)
        print('%s was skipped ...' % enst_id)
        return None
    return transcript_pep


def get_uniprot_aa_seq(uniprot_id):
    """

    Parameters
    ----------
    uniprot_id

    Returns
    -------

    """
    from urllib.request import HTTPError
    fasta_url = 'https://www.uniprot.org/uniprot/' + uniprot_id + '.fasta'
    try:
        url_stream = urllib.request.urlopen(fasta_url)
    except HTTPError:
        return None
    aa_seq = ''
    for l in url_stream.read().decode('ascii').split('\n'):
        if not l.startswith('>'):
            aa_seq += l.strip()
    return aa_seq


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
    dict
        A dictionary where the key is amino acid position and the value is
        the number of variants at this position. One dictionary for missense
        variants and one dictionary for synonymous variants.

    """
    #
    missense_counts = defaultdict(int)
    synonymous_counts = defaultdict(int)
    for variant in variants:
        w = variant[0]  # wild-type amino acid
        v = variant[-1]  # mutant amino acid
        pos = variant[1:-1]  # position in the protein sequence
        if w != v:  # missense variant
            missense_counts[int(pos)] += 1
        else:  # synonymous variant
            synonymous_counts[int(pos)] += 1
    return missense_counts, synonymous_counts


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
    with gzip.open(configs['ensembl_cds'], 'rt') as cds_handle:
        ensembl_cds_dict = SeqIO.to_dict(SeqIO.parse(cds_handle, format='fasta'),
                                 key_function=get_ensembl_accession)

    # CCDS concensus coding sequences
    with gzip.open(configs['ccds_cds'], 'rt') as ccds_handle:
        ccds_dict = SeqIO.to_dict(SeqIO.parse(ccds_handle, format='fasta'),
                                 key_function=get_ccds_accession)
    # ENSEMBL peptide sequences
    with gzip.open(configs['ensembl_pep'], 'rt') as pep_handle:
        pep_dict = SeqIO.to_dict(SeqIO.parse(pep_handle, format='fasta'),
                                 key_function=get_ensembl_accession)

    # parse gnomad transcript-level variants
    with open(configs['gnomad_variants'], 'rt') as variant_handle:
        # transcript_variants will be a dict of dicts where major version
        # ENSEMBL transcript IDs are the first level keys and "ccds", "ensp",
        # "swissprot", "variants" are the second level keys. The value of each
        # second-level key is a Python list.
        transcript_variants = json.load(variant_handle)

    # parse UniProt to PDB chain mapping
    # with open(configs['uniprot_to_pdb'], 'rt') as ipf:
    #    uniprot_to_pdb = dict([l.strip().split() for l in ipf])

    with open(configs['enst_to_pdb'], 'rt') as ipf:
        enst_to_pdb = json.load(ipf)

    output_dir = os.path.abspath(configs['output_dir'])

    # create SIFTS mapping table
    sifts_residue_mapping = SIFTS(configs['sifts_residue_mapping_file'])

    if args.mtr1d:
        suffix = '_mtr1ds.tsv'
    else:
        suffix = '_mtr3ds.tsv'

    # compute the MRT scores for each transcript
    with open(args.transcripts, 'rt') as ipf:
        for transcript in ipf:
            transcript = transcript.strip()
            print('Processing transcript %s' % transcript)
            mtr3ds = []
            contact_stats = []
            if os.path.exists(os.path.join(output_dir, transcript + suffix)) and not args.overwrite:
                print('Scores for %s already exist, skipped.' % transcript)
                continue
            # get the amino acid sequence of the transcript
            try:
                ensp_id = transcript_variants[transcript]['ensp'][0]
            except KeyError:
                logging.critical(
                    'Transcript %s not found in %s',
                    transcript, configs['gnomad_variants']
            )
                continue

            transcript_pep = get_transcript_pep_seq(
                transcript, ensp_id, pep_dict
            )

            # uniprot id for retrieving the corresponding PDB chain
            uniprot_ids = transcript_variants[transcript]['swissprot']
            if len(uniprot_ids) > 1:
                logging.critical(
                    'More than one UniProt IDs found for transcript %s: %s',
                    transcript, ','.join(uniprot_ids)
                )
                sys.exit(1)
            uniprot_id = uniprot_ids[0]

            if transcript_pep is None:
                print('No peptide sequence found for %s' % transcript)

                # fall back to UniProt sequence
                uniprot_pep = get_uniprot_aa_seq(uniprot_id)
                if uniprot_pep is None:
                    continue
                else:
                    transcript_pep = uniprot_pep

            # get all variants of this transcript reported in gnomAD
            try:
                variants = transcript_variants[transcript]['variants']
            except KeyError:
                logging.critical('No variants found in %s in gnomAD', transcript)
                logging.critical('%s was skipped ...', transcript)
                continue

            # tabulate variants at each site
            missense_counts, synonymous_counts = count_variants(variants)

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

            # skip if the CDS is incomplete
            if len(transcript_cds) / 3 != len(transcript_pep) + 1:
                print('Incomplete CDS for', transcript, '. Skipped.')
                continue

            # check that the CDS does not contain invalid nucleotides
            if any([x not in {'A', 'T', 'C', 'G'} for x in set(transcript_cds)]):
                print('Invalid CDS! Skipped.')
                print(transcript_cds)
                continue

            # calculate expected counts for each codon
            expected_counts = count_cds_ns(transcript_cds)

            # only compute MTR1D scores if asked on the command-line
            if args.mtr1d:
                for i, a in enumerate(transcript_pep, start=1):
                    mtr3ds.append([transcript, ensp_id, uniprot_id, i, a] +
                        list(compute_mtr1d(i, missense_counts, synonymous_counts,
                                      expected_counts, window=31))
                        )
                with open(file=os.path.join(output_dir, transcript + suffix), mode='wt') as opf:
                    header = ['TRANSCRIPT', 'PROTEIN', 'UNIPROT', 'ENSP_POS', 'ENSP_AA',
                    'MISSENSE', 'SYNONYMOUS', 'MTR1D']
                    csv_writer = csv.writer(opf, delimiter='\t')
                    csv_writer.writerow(header)
                    csv_writer.writerows(mtr3ds)
                print('Finished calcualting MTR1D scores for %s!' % transcript)
                continue

            # get the PDB ID and PDB chain associated with this transcript
            try:
                pdb_id, pdb_chain = enst_to_pdb[transcript]
            except KeyError:
                logging.critical(
                    '%s not found in given database %s',
                    transcript, configs['enst_to_pdb']
                )
                logging.critical('%s was skipped ...', transcript)
                print('%s not found in given database %s' %
                      (transcript, configs['enst_to_pdb']))
                print('%s was skipped ...' % transcript)
                continue

            # print message
            print('Estimating MRT3D scores for:', transcript, ensp_id,
                  uniprot_id, pdb_id, pdb_chain)

            chain = pdb_utils.get_pdb_chain(
                pdb_id, pdb_chain, configs['pdb_dir']
            )

            if chain is None:
                print('No chain was available for %s' % transcript)
                continue

            all_aa_residues = [aa for aa in chain.get_residues() if is_aa(aa)]
            all_contacts = search_for_all_contacts(all_aa_residues, radius=4)

            # index all contacts by residue ID
            indexed_contacts = defaultdict(list)
            for c in all_contacts:
                indexed_contacts[c.get_res_a()].append(c.get_res_b())
                indexed_contacts[c.get_res_b()].append(c.get_res_a())


            uniprot_to_pdb_mapping = sifts_residue_mapping.uniprot_to_pdb(
                uniprot_id, pdb_id, pdb_chain
            )
            pdb_to_uniprot_mapping = sifts_residue_mapping.pdb_to_uniprot(
                pdb_id, pdb_chain, uniprot_id
            )

            # 
            id_fields = [transcript, ensp_id, uniprot_id, pdb_id, pdb_chain]
            
            for i, a in enumerate(transcript_pep, start=1):
                # @TODO need to get the corresponding position in the PDB file
                # It is highly likely that the amino acid at position i in the
                # PDB file will not be the same amino acid at position i in the
                # ENSP peptide sequence
                try:
                    pdb_pos = uniprot_to_pdb_mapping[i]
                    res = chain[pdb_pos]
                except KeyError:
                    # mtr3d.append((i, compute_mtr1d(i, missense_counts, synonymous_counts,
                    #                            expected_counts, window=31)))
                    # print('Residue', i, 'not found in', pdb_file)
                    logging.critical(
                        'Residue %s in %s not found in chain %s in PDB file: %s', 
                        i, ensp_id, pdb_chain, pdb_id
                    )
                    mtr3ds.append(id_fields + [i, a] + [np.nan] * 7)
                    # contact_stats.append(
                    #    (transcript, ensp_id, pdb_id, pdb_chain, i, np.nan)
                    # )
                    continue
                    
                # check that the amino acid in ENSP sequence matches 
                # that in the PDB structure
                pdb_aa = seq1(res.get_resname())
                if a != pdb_aa:
                    logging.critical(
                        'Residue in %s and %s of %s did not match that in PDB %s '
                        'chain %s at %s: %s vs %s',
                        uniprot_id, ensp_id, transcript, pdb_id, pdb_chain, i, a, pdb_aa
                    )
                    logging.critical(
                        'Please first check UniProt sequence is identical to ENSP '
                        'sequence. If this is true, check if there is any oddity '
                        'in the SIFTS residue-level mapping.'
                    )
                    mtr3ds.append(id_fields + [i, a, pdb_pos, pdb_aa] + [np.nan] * 5)
                    # contact_stats.append(
                    #    (transcript, ensp_id, pdb_id, pdb_chain, i, np.nan)
                    # )
                    continue
                    # sys.exit(1)

                contact_res = indexed_contacts[res]

                # exclude contacting residues that are adjacent in sequence
                # contacts_pdb_pos = [
                #     r.get_id()[1] for r in contact_res
                #     if r.get_id()[1] - i > 3 or r.get_id()[1] - i < -3
                # ]
                contacts_pdb_pos = [r.get_id()[1] for r in contact_res]
                
                seq_seps = ';'.join(str(x) for x in 
                                    [i - pdb_pos for i in contacts_pdb_pos])

                # collect statistics about the count of contacts
                # contact_stats.append(
                #    (transcript, ensp_id, pdb_id, pdb_chain, i, len(contacts_pdb_pos))
                # )
                num_contacts = len(contacts_pdb_pos)

                total_missense_obs = missense_counts.setdefault(i, 0)
                total_synonymous_obs = synonymous_counts.setdefault(i, 0)
                try:
                    total_missense_exp = expected_counts[i - 1][0]
                    total_synonymous_exp = expected_counts[i - 1][1]
                except IndexError:
                    print('list index out of range:', i)
                if contacts_pdb_pos:
                    for j in contacts_pdb_pos:
                        # @TODO need to get the corresponding position in ENSP
                        # amino acid sequence using the SIFTS mapping table
                        try:
                            ensp_pos = pdb_to_uniprot_mapping[j]
                        except KeyError:
                            logging.critical(
                                'PDB residue %s in %s chain %s not found '
                                'in %s, skipped',
                                j, pdb_id, pdb_chain, ensp_id
                            )
                            continue

                        # count the total # observed variants in contacting residues
                        total_missense_obs += missense_counts.setdefault(ensp_pos, 0)
                        total_synonymous_obs += synonymous_counts.setdefault(ensp_pos, 0)
                        # count the total # expected variants
                        try:
                            total_missense_exp += expected_counts[ensp_pos - 1][0]
                            total_synonymous_exp += expected_counts[ensp_pos - 1][1]
                        except IndexError:
                            logging.critical(
                                'PDB residue %s in %s chain %s out of the range '
                                'of %s codons',
                                j, pdb_id, pdb_chain, transcript
                            )
                            break
                    # compute the fraction of observed missense variants
                    try:
                        missense_frac_obs = total_missense_obs / (total_missense_obs + total_synonymous_obs)
                    except ZeroDivisionError:
                        logging.critical(
                            'No variants were observed in the neighborhood of '
                            'position %s in %s of %s in gnomAD', i, ensp_id, transcript)
                        logging.critical(
                            'The MTR3D score for this position was set to 0.'
                        )
                        mtr3ds.append(id_fields + 
                            [i, a, pdb_pos, pdb_aa] + 
                            [seq_seps, num_contacts, total_missense_obs, 
                             total_synonymous_obs, 0.0]
                        )
                        continue

                    # compute the fraction of expected missense variants
                    missense_frac_exp = total_missense_exp / (total_missense_exp + total_synonymous_exp)
                    mtr3ds.append(id_fields + [i, a, pdb_pos, pdb_aa] + 
                        [seq_seps, num_contacts, total_missense_obs, 
                         total_synonymous_obs, missense_frac_obs / missense_frac_exp])
                else:
                    logging.critical(
                        'No contacts were found in the neighborhood of '
                        'position %s in the PDB file %s chain %s of mapped '
                        'residue %s in %s of %s',
                        pdb_pos, pdb_id, pdb_chain, i, ensp_id, transcript
                    )
                    logging.critical(
                        'The MTR score for this position was set to 1.0'
                    )
                    mtr3ds.append(id_fields + [i, a, pdb_pos, pdb_aa] + 
                                  [seq_seps, num_contacts, 0, 0, 1.0])

            with open(file=os.path.join(output_dir, transcript + '_mtr3ds.tsv'), mode='wt') as opf:
                header = ['TRANSCRIPT', 'PROTEIN', 'UNIPROT', 'PDB', 'CHAIN',
                    'ENSP_POS', 'ENSP_AA', 'PDB_POS', 'PDB_AA', 'CONTACS',
                    'MISSENSE', 'SYNONYMOUS', 'MTR3D']
                csv_writer = csv.writer(opf, delimiter='\t')
                csv_writer.writerow(header)
                csv_writer.writerows(mtr3ds)

            # with open(file=transcript + '_contact_stats.tsv', mode='wt') as opf:
            #    csv_writer = csv.writer(opf, delimiter='\t')
            #    csv_writer.writerows(contact_stats)


if __name__ == '__main__':
    main()
