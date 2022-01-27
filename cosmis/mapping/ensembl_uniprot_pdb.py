#!/usr/bin/env python3

import os
import pandas as pd
import urllib
from cosmis.mapping.sifts import SIFTS
from cosmis.utils import pdb_utils
from collections import defaultdict


SIFTS_ENSEMBL_URL = 'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/' \
                    'tsv/pdb_chain_ensembl.tsv.gz'


class EnsemblUniProtPDB:
    """

    """
    def __init__(self, sifts_ensembl=None, sifts_uniprot=None, pdb_path=None):
        """

        Parameters
        ----------
        sifts_uniprot : str
            Residue-level mapping file.

        sifts_ensembl : str
            Chain-level mapping file.
        """
        if sifts_ensembl is None:
            local_path = os.path.abspath(
                './pdb_chain_ensembl.tsv.gz'
            )
            if os.path.exists(local_path):
                sifts_ensembl = local_path
            else:
                # download sifts mapping file
                urllib.request.urlretrieve(SIFTS_ENSEMBL_URL, local_path)
                sifts_ensembl = local_path

        self.mapping_table = self._create_mapping_table(sifts_ensembl)
        self.sifts_uniprot = sifts_uniprot
        self.pdb_path = pdb_path

    def _create_mapping_table(self, sifts_mapping_file):
        """

        Parameters
        ----------
        sifts_mapping_file

        Returns
        -------

        """
        # load the mapping file into a Pandas DataFrame
        sifts_table = pd.read_csv(
            sifts_mapping_file,
            sep='\t',
            compression='gzip',
            comment='#',
            na_values='None',
            low_memory=False
        )

        sifts_table.rename(
            inplace=True,
            columns={
                'PDB': 'pdb_id',
                'CHAIN': 'pdb_chain',
                'SP_PRIMARY': 'uniprot_id',
                'GENE_ID': 'ensg_id',
                'TRANSCRIPT_ID': 'enst_id',
                'TRANSLATION_ID': 'ensp_id',
                'EXON_ID': 'exon_id'
            }
        )

        return sifts_table

    def enst_to_pdb(self, enst_id, uniprot_id=None):
        """

        Parameters
        ----------
        enst_id
        uniprot_id

        Returns
        -------

        """
        enst_id = enst_id.upper()
        query_str = 'enst_id == ' + '"' + enst_id + '"'

        if uniprot_id is not None:
            query_str += ' and uniprot_id == ' + '"' + uniprot_id + '"'

        # query the mapping table
        hits = self.mapping_table.query(query_str)

        if hits.empty:
            return None, None

        pdb_chains = defaultdict(list)
        for _, r in hits.iterrows():
            # skip records where PDB ID or CHAIN ID is empty string
            if r['pdb_id'] and r['pdb_chain']:
                pdb_chains[r['pdb_id']].append(r['pdb_chain'])

        # remove duplicates but keep ordering
        uniq_pdb_chains = []
        for k, v in pdb_chains.items():
            uniq_pdb_chains.append((k, v[0]))

        if not uniq_pdb_chains:
            return None, None

        # return the pdb chain that has the largest coverage
        # of the transcript protein sequence, this could be implemented
        # ad hoc, but for now, we rely on SIFTS residue-level mapping
        sifts = SIFTS(sifts_uniprot=self.sifts_uniprot, xml_dir=self.pdb_path)
        max_len = 0
        best_resolution = pdb_utils.get_resolution(uniq_pdb_chains[0][0], self.pdb_path)
        best_pdb_id = ''
        best_chain_id = ''
        for pdb_id, chain_id in uniq_pdb_chains:
            if not (pdb_id and chain_id):
                print('PDB ID is empty string for', enst_id, ', skipped')
                continue
            residue_mapping = sifts.pdb_to_uniprot(pdb_id, chain_id)
            if residue_mapping is None:
                print('Failed to obtained residue mapping from SIFTS xml file.')
                continue
            # check sequence coverage
            if max_len < len(residue_mapping):
                max_len = len(residue_mapping)
                best_pdb_id = pdb_id
                best_chain_id = chain_id
            # check resolution
            elif max_len == len(residue_mapping):
                resolution = pdb_utils.get_resolution(pdb_id, self.pdb_path)
                if best_resolution is None:
                    best_resolution = resolution
                if resolution is not None and resolution < best_resolution:
                    best_pdb_id = pdb_id
                    best_chain_id = chain_id
                    best_resolution = resolution

        return best_pdb_id, best_chain_id

    def enst_to_uniprot(self, enst_id):
        """

        Parameters
        ----------
        enst_id

        Returns
        -------

        """
        enst_id = enst_id.upper()
        query_str = 'enst_id == ' + '"' + enst_id + '"'
        hits = self.mapping_table.query(query_str)
        return set(hits.unique())

    def uniprot_to_pdb(self, uniprot_id, multimeric_state=0):
        """

        Parameters
        ----------
        uniprot_id

        Returns
        -------

        """
        uniprot_id = uniprot_id.upper()
        query_str = 'uniprot_id == ' + '"' + uniprot_id + '"'

        # query the mapping table
        hits = self.mapping_table.query(query_str)

        if hits.empty:
            return None, None

        pdb_chains = defaultdict(set)
        for _, r in hits.iterrows():
            # skip records where PDB ID or CHAIN ID is empty string
            if r['pdb_id'] and r['pdb_chain']:
                pdb_chains[r['pdb_id']].add(r['pdb_chain'])

        uniq_pdb_chains = []
        for k, v in pdb_chains.items():
            if multimeric_state:
                if len(v) != multimeric_state:
                    continue

                # make sure that the PDB file has exactly multimeric_state
                # number of chains
                query_by_pdb = 'pdb_id == ' + '"' + k + '"'
                pdb_hits = self.mapping_table.query(query_by_pdb)
                chain_records = set()
                for _, x in pdb_hits.iterrows():
                    if x['pdb_id'] and x['pdb_chain']:
                        chain_records.add(x['pdb_id'] + x['pdb_chain'])
                if len(chain_records) != multimeric_state:
                    continue

            # remove duplicates but keep ordering
            uniq_pdb_chains.append((k, list(v)[0]))

        if not uniq_pdb_chains:
            return None, None

        # return the pdb chain that has the largest coverage
        # of the transcript protein sequence, this could be implemented
        # ad hoc, but for now, we rely on SIFTS residue-level mapping
        sifts = SIFTS(sifts_uniprot=self.sifts_uniprot, xml_dir=self.pdb_path)
        max_len = 0
        best_resolution = 5.0
        best_pdb_id = ''
        best_chain_id = ''
        for pdb_id, chain_id in uniq_pdb_chains:
            if not (pdb_id and chain_id):
                print('PDB ID is empty string for', uniprot_id, ', skipped')
                continue
            # skip if resolution worse than 5 angstrom
            resolution = pdb_utils.get_resolution(pdb_id, self.pdb_path)
            # if resolution is None, then it is likely an NMR model
            # if resolution is None:
            #     print('{} {} is an NMR model. EXCLUDED!'.format(pdb_id, chain_id))
            #    continue
            if resolution is not None and resolution > 5.0:
                print('{} {} resolution > 5 angstrom. EXCLUDED!'.format(pdb_id, chain_id))
                continue
            residue_mapping = sifts.pdb_to_uniprot(pdb_id, chain_id, uniprot_id)
            if residue_mapping is None:
                print('Failed to obtain residue mapping from SIFTS xml file.')
                continue
            # check sequence coverage
            if max_len < len(residue_mapping):
                max_len = len(residue_mapping)
                best_pdb_id = pdb_id
                best_chain_id = chain_id
                if resolution is not None:
                    best_resolution = resolution
            # check resolution
            elif max_len == len(residue_mapping):
                if resolution is not None and resolution < best_resolution:
                    best_pdb_id = pdb_id
                    best_chain_id = chain_id
                    best_resolution = resolution

        return best_pdb_id, best_chain_id

    def pdb_to_uniprot(self, pdb_id, multimer=False):
        """

        Parameters
        ----------
        pdb_id
        multimer

        Returns
        -------

        """
        pass


def main():
    """

    Returns
    -------

    """
    ensembl_uniprot_pdb = EnsemblUniProtPDB()

    test_uniprot_id = 'Q9Y4Z0'

    pdb_id, chain_id = ensembl_uniprot_pdb.uniprot_to_pdb(test_uniprot_id)

    print(pdb_id, chain_id)


if __name__ == '__main__':
    main()
