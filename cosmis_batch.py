#!/usr/bin/env python3

import os, csv
import gzip
import json
import numpy as np
from collections import defaultdict
from argparse import ArgumentParser
from cosmis.utils import pdb_utils, seq_utils
from Bio import SeqIO
from Bio.SeqUtils import seq1
from Bio.PDB import PDBParser, is_aa


def parse_cmd():
    """

    Returns
    -------

    """
    parser = ArgumentParser()
    parser.add_argument('-c', '--config', dest='config', required=True,
                        type=str, help='A JSON file specifying options.')
    parser.add_argument('-i', '--input', dest='input', required=True,
                        type=str, help='''Ensembl transcript ID of the 
                        transcript for which to compute a MTR3D profile.''')
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
    try:
        transcript_pep = pep_dict[ensp_id].seq
    except KeyError:
        print('%s not found in given database' % ensp_id)
        print('%s was skipped ...' % enst_id)
        return None
    return transcript_pep


def get_pdb_chain(pdb_file, pdb_chain):
    """
    Creates a Bio.PDB.Chain object for the requested PDB chain.

    Parameters
    ----------
    pdb_id : str
        Identifier of the PDB chain as a five-letter string.
    pdb_chain : str
        Identifier of the PDB chain as a five-letter string.
    pdb_db : str
        Path to the local PDB database.

    Returns
    -------
    Bio.PDB.Chain
        The requested PDB chain as a Bio.PDB.Chain object.

    """
    # read in the PDB file
    pdb_parser = PDBParser(PERMISSIVE=1)
    try:
        structure = pdb_parser.get_structure(id="given_pdb", file=pdb_file)
    except (FileNotFoundError, ValueError) as e:
        print('PDB file cannot be retrieved', pdb_file)
        return None
    try:
        chain = structure[0][pdb_chain]
    except KeyError:
        print('No chain ' + pdb_chain + ' was found in ' + pdb_file)
        return None
    return chain


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


def main():
    """

    Returns
    -------

    """
    # parse command-line arguments
    args = parse_cmd()

    # parse configuration file
    configs = parse_config(args.config)

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

    # ENSEMBL peptide sequences
    print('Reading ENSEMBL protein sequence database ...')
    with gzip.open(configs['ensembl_pep'], 'rt') as pep_handle:
        pep_dict = SeqIO.to_dict(
            SeqIO.parse(pep_handle, format='fasta'),
            key_function=get_ensembl_accession
        )

    # parse gnomad transcript-level variants
    print('Reading gnomAD variant database ...')
    with open(configs['gnomad_variants'], 'rt') as variant_handle:
        # transcript_variants will be a dict of dicts where major version
        # ENSEMBL transcript IDs are the first level keys and "ccds", "ensp",
        # "swissprot", "variants" are the second level keys. The value of each
        # second-level key is a Python list.
        transcript_variants = json.load(variant_handle)

    # get phylop scores
    print('Loading PhyloP scores ...')
    with gzip.open(configs['enst_to_phylop'], 'rt') as ipf:
        enst_to_phylop = json.load(ipf)

    # directory where to store the output files
    output_dir = os.path.abspath(configs['output_dir'])

    # read transcript to swiss model mapping
    with open(args.input, 'rt') as ipf:
        transcripts = [line.strip().split() for line in ipf]

    # compute the MRT scores
    for x, y in transcripts:
        transcript = x
        pdb_file = 'data/swiss-model/SWISS-MODEL_Repository/' + y + '.pdb'
        pdb_chain = y[-1]

        if os.path.exists(
                os.path.join(output_dir, transcript + '_cosmis.tsv')
        ) and not args.overwrite:
            continue

        print('Processing transcript %s' % transcript)

        cosmis = []

        # get the amino acid sequence of the transcript
        try:
            # Ensembl peptide ID for the transcript
            ensp_id = transcript_variants[transcript]['ensp'][0]
        except KeyError:
            print(
                'Transcript %s not found in %s' %
                (transcript, configs['gnomad_variants'])
            )
            print('Skip to the next transcript...')
            continue

        # get the peptide sequence from peptide sequence database
        transcript_pep_seq = get_transcript_pep_seq(
            transcript, ensp_id, pep_dict
        )

        # get all variants of this transcript reported in gnomAD
        try:
            variants = transcript_variants[transcript]['variants']
        except KeyError:
            print('No variants found for %s in gnomAD', transcript)
            print('Skip to the next transcript...')
            continue

        # get the coding sequence of the transcript
        try:
            transcript_cds = ensembl_cds_dict[transcript].seq
        except KeyError:
            print(
                '''No CDS found in Ensembl CDS database! 
                Looking for it in the CCDS database ...'''
            )
            transcript_cds = None

        if transcript_cds is None:
            try:
                ccds_id = transcript_variants[transcript]['ccds'][0]
                transcript_cds = ccds_dict[ccds_id].seq
            except KeyError:
                print('ERROR: No CDS found in CCDS database!')
                print('Skip to the next transcript...')
                continue

        if transcript_pep_seq is None:
            print('ERROR: No peptide sequence found in Ensembl database!')
            print('Skip to the next transcript...')
            continue

        # stop if the CDS is invalid
        if len(transcript_cds) / 3 != len(transcript_pep_seq) + 1:
            print('ERROR: incomplete CDS for', transcript)
            print('Skip to the next transcript...')
            continue

        # check that the CDS does not contain invalid nucleotides
        if any([x not in {'A', 'T', 'C', 'G'} for x in set(transcript_cds)]):
            print('ERROR: invalid CDS for', transcript_cds)
            print('Skip to the next transcript...')
            continue

        # get the phyloP scores for the current transcript
        try:
            transcript_phylop_scores = enst_to_phylop[transcript]['phylop']
        except KeyError:
            print('No phyloP scores are available for {}'.format(transcript))
            continue

        # print message
        print('Computing COSMIS features for:', transcript, ensp_id,pdb_file)

        chain = get_pdb_chain(pdb_file, pdb_chain)

        if chain is None:
            print(
                'ERROR: %s not found in structure: %s!' % (pdb_chain, pdb_file))
            print('Skip to the next transcript...')
            continue

        all_aa_residues = [aa for aa in chain.get_residues() if is_aa(aa)]
        all_contacts = pdb_utils.search_for_all_contacts(all_aa_residues,
                                                         radius=8)

        # calculate expected counts for each codon
        codon_mutation_rates = seq_utils.get_codon_mutation_rates(
            transcript_cds)
        all_cds_ns_counts = seq_utils.count_poss_ns_variants(transcript_cds)
        cds_ns_sites = seq_utils.count_ns_sites(transcript_cds)

        if len(codon_mutation_rates) < len(all_aa_residues):
            print('ERROR: peptide sequence has less residues than structure!')
            print('Skip to the next transcript...')
            continue

        # tabulate variants at each site
        # missense_counts and synonymous_counts are dictionary that maps
        # amino acid positions to variant counts
        missense_counts, synonymous_counts = count_variants(variants)

        # convert variant count to site variability
        site_variability_missense = {
            pos: 1 for pos, _ in missense_counts.items()
        }
        site_variability_synonymous = {
            pos: 1 for pos, _ in synonymous_counts.items()
        }

        # compute the total number of missense variants
        total_mis_counts = 0
        for k, v in missense_counts.items():
            total_mis_counts += v

        # permutation test
        permuted_missense_mutations = seq_utils.permute_missense(
            total_mis_counts, len(transcript_pep_seq)
        )

        # index all contacts by residue ID
        indexed_contacts = defaultdict(list)
        for c in all_contacts:
            indexed_contacts[c.get_res_a()].append(c.get_res_b())
            indexed_contacts[c.get_res_b()].append(c.get_res_a())

        valid_case = True
        for seq_pos, seq_aa in enumerate(transcript_pep_seq, start=1):
            # check that the amino acid in ENSP sequence matches 
            # that in the PDB structure
            try:
                res = chain[seq_pos]
            except KeyError:
                print('PDB file is missing residue:', seq_aa, 'at', seq_pos)
                continue
            pdb_aa = seq1(res.get_resname())
            if seq_aa != pdb_aa:
                print('Residue in ENSP did not match that in PDB at', seq_pos)
                print('Skip to the next transcript...')
                valid_case = False
                break

            contact_res = indexed_contacts[res]
            num_contacts = len(contact_res)
            contacts_pdb_pos = [r.get_id()[1] for r in contact_res]
            seq_seps = ';'.join(
                str(x) for x in [i - seq_pos for i in contacts_pdb_pos]
            )

            mis_var_sites = site_variability_missense.setdefault(seq_pos, 0)
            total_mis_sites = cds_ns_sites[seq_pos - 1][0]
            syn_var_sites = site_variability_synonymous.setdefault(seq_pos, 0)
            total_syn_sites = cds_ns_sites[seq_pos - 1][1]
            total_missense_obs = missense_counts.setdefault(seq_pos, 0)
            total_synonymous_obs = synonymous_counts.setdefault(seq_pos, 0)
            total_missense_poss = all_cds_ns_counts[seq_pos - 1][0]
            total_synonyms_poss = all_cds_ns_counts[seq_pos - 1][1]
            total_synonymous_rate = codon_mutation_rates[seq_pos - 1][0]
            total_missense_rate = codon_mutation_rates[seq_pos - 1][1]
            phylop_scores = transcript_phylop_scores[seq_pos - 1][3]
            for j in contacts_pdb_pos:
                # count the total # observed variants in contacting residues
                mis_var_sites += site_variability_missense.setdefault(j, 0)
                syn_var_sites += site_variability_synonymous.setdefault(j, 0)
                total_missense_obs += missense_counts.setdefault(j, 0)
                total_synonymous_obs += synonymous_counts.setdefault(j, 0)

                # count the total # expected variants
                try:
                    total_missense_poss += all_cds_ns_counts[j - 1][0]
                    total_synonyms_poss += all_cds_ns_counts[j - 1][1]
                    total_synonymous_rate += codon_mutation_rates[j - 1][0]
                    total_missense_rate += codon_mutation_rates[j - 1][1]
                    total_mis_sites += cds_ns_sites[j - 1][0]
                    total_syn_sites += cds_ns_sites[j - 1][1]
                except IndexError:
                    valid_case = False
                    break
            if not valid_case:
                break

            try:
                seq_context = seq_utils.get_codon_seq_context(
                    contacts_pdb_pos + [seq_pos], transcript_cds
                )
            except IndexError:
                break

            # compute the GC content of the sequence context
            if len(seq_context) == 0:
                print('No nucleotides were found in sequence context!')
                continue
            gc_fraction = seq_utils.gc_content(seq_context)

            contact_res_indices = [pos - 1 for pos in contacts_pdb_pos] + [
                seq_pos - 1]
            n = np.sum(permuted_missense_mutations[:, contact_res_indices].sum(axis=1) <= total_mis_sites)
            p_value = n / 10000

            # compute the fraction of expected missense variants
            cosmis.append(
                [
                    transcript, ensp_id, seq_pos, seq_aa, seq_seps,
                    num_contacts + 1,
                    syn_var_sites,
                    '{:.3f}'.format(total_syn_sites),
                    mis_var_sites,
                    '{:.3f}'.format(total_mis_sites),
                    total_synonyms_poss,
                    total_missense_poss,
                    '{:.3f}'.format(gc_fraction),
                    '{:.3e}'.format(total_synonymous_rate),
                    total_synonymous_obs,
                    '{:.3e}'.format(total_missense_rate),
                    total_missense_obs,
                    p_value,
                    '{:.3f}'.format(np.mean(phylop_scores)),
                    total_mis_counts,
                    len(transcript_pep_seq)
                ]
            )

        if not valid_case:
            continue

        with open(
            file=os.path.join(output_dir, transcript + '_cosmis.tsv'),
            mode='wt'
        ) as opf:
            header = [
                'enst_id', 'ensp_id', 'ensp_pos', 'ensp_aa',
                'seq_separations', 'num_contacts',
                'syn_var_sites', 'total_syn_sites',
                'mis_var_sites', 'total_mis_sites',
                'synonymous_poss', 'missense_poss', 'gc_content',
                'synonymous_rate', 'synonymous_obs', 'missense_rate',
                'missense_obs', 'p_value',
                'phylop_score', 'enst_mis_counts', 'ensp_length'
            ]
            csv_writer = csv.writer(opf, delimiter='\t')
            csv_writer.writerow(header)
            csv_writer.writerows(cosmis)


if __name__ == '__main__':
    main()
