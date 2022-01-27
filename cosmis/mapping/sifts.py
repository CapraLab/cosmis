#!/usr/bin/env python3

import os
import gzip
import pandas as pd
import numpy as np
import urllib
import wget
import signal
import xml.etree.ElementTree as ET


SIFTS_URL = 'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/' \
            'pdb_chain_uniprot.tsv.gz'
SIFTS_XML_SCHEMA = 'http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd'


def wget_timeout_handler(signum, frame):
    """
    

    Parameters
    ----------
    signum : TYPE
        DESCRIPTION.
    frame : TYPE
        DESCRIPTION.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    None.

    """
    print('WGET TIMEOUT!')
    raise Exception('wget taking too long!')    


class XMLNamespaces:
    """
    A utility class to enable cleaner parse of XML documents with namespaces.
    """
    def __init__(self, **kwargs):
        """

        Parameters
        ----------
        kwargs
        """
        self.namespaces = {}
        for name, uri in kwargs.items():
            self.register(name, uri)

    def register(self, name, uri):
        """

        Parameters
        ----------
        name
        uri

        Returns
        -------

        """
        self.namespaces[name] = '{' + uri + '}'

    def __call__(self, path):
        """

        Parameters
        ----------
        path

        Returns
        -------

        """
        return  path.format_map(self.namespaces)


class SIFTS:
    """
    Provide position mapping between PDB and UniProt based on SIFTS
    mapping table.
    """
    def __init__(self, sifts_uniprot=None, xml_dir=None):
        """

        Parameters
        ----------
        sifts_uniprot : str
        """
        if sifts_uniprot is None:
            print(
                'SIFTS mapping file not given, now downloading it from',
                SIFTS_URL
            )
            # first try to look for the file locally
            local_path = os.path.abspath(
                './pdb_chain_uniprot.tsv.gz'
            )
            if os.path.exists(local_path):
                sifts_uniprot = local_path
            else:
                # download the mapping file to local path
                urllib.request.urlretrieve(SIFTS_URL, local_path)
                sifts_uniprot = local_path
        
        self.sifts_table = self._create_mapping_table(sifts_uniprot)
        
        if xml_dir is None:
            self.xml_dir = os.path.abspath('/tmp/')
        else:
            self.xml_dir = xml_dir

    def _create_mapping_table(self, sifts_uniprot):
        """

        Parameters
        ----------
        sifts_uniprot : str
            Path to the computed mapping table.

        Returns
        -------

        """
        sifts_table = pd.read_csv(
            sifts_uniprot,
            sep='\t',
            comment='#',
            compression='gzip',
            na_values='None',
            low_memory=False
        )

        sifts_table.rename(
            inplace=True,
            columns={
                'PDB': 'pdb_id',
                'CHAIN': 'pdb_chain',
                'SP_PRIMARY': 'uniprot_id',
                'RES_BEG': 'resseq_beg',
                'RES_END': 'resseq_end',
                'PDB_BEG': 'pdb_beg',
                'PDB_END': 'pdb_end',
                'SP_BEG': 'uniprot_beg',
                'SP_END': 'uniprot_end'
            }
        )

        return sifts_table

    def pdb_to_uniprot(self, pdb_id, pdb_chain=None, uniprot_id=None):
        """

        Parameters
        ----------
        pdb_id : str
            4-letter PDB identifier
        pdb_chain : str, optional (default: None)
            PDB chain identifier (if not given, all
            chains for the PDB entry will be returned).
        uniprot_id : str, optional (default: None)

        Returns
        -------

        """
        pdb_id = pdb_id.lower()
        query_str = 'pdb_id == ' + '"' + pdb_id + '"'

        # filter by PDB chain
        if pdb_chain is not None:
            query_str += ' and pdb_chain == ' + '"' + pdb_chain + '"'

        # filter by UniProt ID
        if uniprot_id is not None:
            query_str += ' and uniprot_id == ' + '"' + uniprot_id + '"'

        hits = self.sifts_table.query(query_str)

        # create position mapping from PDB to UniProt
        mapping = {}
        for _, r in hits.iterrows():
            # try to cast PDB positions to Python int objects
            if str(r['pdb_beg']).isdigit():
                pdb_beg = int(r['pdb_beg'])
            else:
                pdb_beg = np.nan
            if str(r['pdb_end']).isdigit():
                pdb_end = int(r['pdb_end'])
            else:
                pdb_end = np.nan

            uniprot_beg, uniprot_end = r['uniprot_beg'], r['uniprot_end']

            if (not pd.isna(pdb_beg)) and (not pd.isna(pdb_end)) and \
            pdb_end - pdb_beg == uniprot_end - uniprot_beg:
                mapping.update(
                    {x: y for x, y in zip(
                        range(pdb_beg, pdb_end + 1),
                        range(uniprot_beg, uniprot_end + 1)
                    )}
                )
            else:  # create mapping from SIFTS XML mapping file
                mapping = self.pdb_to_uniprot_xml(pdb_id, pdb_chain, uniprot_id)
                break

        return mapping

    def uniprot_to_pdb(self, uniprot_id, pdb_id=None, pdb_chain=None):
        """

        Parameters
        ----------
        pdb_id
        pdb_chain
        uniprot_id

        Returns
        -------

        """
        # construct the query string
        uniprot_id = uniprot_id.upper()
        query_str = 'uniprot_id == ' + '"' + uniprot_id + '"'

        if pdb_id is not None and pdb_chain is not None:
            reduce_chains = False
            pdb_id = pdb_id.lower()
            query_str += ' and pdb_id == ' + '"' + pdb_id + '"' + \
                         ' and pdb_chain == ' + '"' + pdb_chain + '"'
        else:
            reduce_chains = True

        hits = self.sifts_table.query(query_str)

        # only retain one chain if this option is set to True
        # @TODO this chunk needs to be revisited to make sure that all the
        # segments of the chain are retained
        if reduce_chains:
            hits = hits.groupby('pdb_id').first().reset_index()

        # create position mapping from UniProt to PDB
        mapping = {}
        for _, r in hits.iterrows():
            # try to cast PDB positions to Python int objects
            if str(r['pdb_beg']).isdigit():
                pdb_beg = int(r['pdb_beg'])
            else:
                pdb_beg = np.nan
            if str(r['pdb_end']).isdigit():
                pdb_end = int(r['pdb_end'])
            else:
                pdb_end = np.nan

            uniprot_beg, uniprot_end = r['uniprot_beg'], r['uniprot_end']

            if (not pd.isna(pdb_beg)) and (not pd.isna(pdb_end)) and \
            pdb_end - pdb_beg == uniprot_end - uniprot_beg:
                mapping.update(
                    {y: x for x, y in zip(
                        range(pdb_beg, pdb_end + 1),
                        range(uniprot_beg, uniprot_end + 1)
                    )}
                )
            else:  # create mapping from SIFTS XML mapping file
                pdb_to_uniprot = self.pdb_to_uniprot_xml(pdb_id, pdb_chain, uniprot_id)
                if pdb_to_uniprot is not None:
                    mapping = {
                        v: k for k, v in pdb_to_uniprot.items()
                    }
                else:
                    mapping.update(
                        {x: x for x in range(uniprot_beg, uniprot_end + 1)}
                    )
                break

        return mapping

    def pdb_to_uniprot_xml(
            self, pdb_id, pdb_chain, uniprot_id=None, timeout=600
        ):
        """

        Parameters
        ----------
        pdb_id : str
            Four-letter PDB ID
        pdb_chain : str
            One-letter PDB chain ID
        xml_path : str
            Path to the SIFTS XML residue-level mapping file
        uniprot_id : str
            UniProt ID

        Returns
        -------
        dict
            Mapping of PDB residue numbering to UniProt residue numbering.

        """
        if not os.path.exists(os.path.join(self.xml_dir, pdb_id[1:3])):
            os.mkdir(os.path.join(self.xml_dir, pdb_id[1:3]))
        xml_file = os.path.join(self.xml_dir, pdb_id[1:3], pdb_id + '.xml.gz')
        if not os.path.exists(xml_file):
            xml_file_ftp = 'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/{}.xml.gz'.format(pdb_id) 
            try:
                print('Downloding', xml_file_ftp)
                # timeout if the download is taking too long
                signal.signal(signal.SIGALRM, wget_timeout_handler)
                signal.alarm(timeout)
                wget.download(xml_file_ftp, xml_file)
                # set the alarm off
                signal.alarm(0)
            except:
                print('Cannot retrieve', xml_file)
                return None
        
        try:
            if os.stat(xml_file).st_size == 0:
                print(xml_file, 'is empty!')
                return None
            with gzip.open(xml_file, 'rt') as ipf:
                xml_mapping = ET.parse(ipf)
        except EOFError:
            print('ERROR reading', xml_file)
            return None

        # print('{} already exists locally.'.format(xml_file))
        xml_ns = XMLNamespaces(entry=SIFTS_XML_SCHEMA)
        all_res = list(
            xml_mapping.getroot().iterfind(
                xml_ns('{entry}entity/{entry}segment/{entry}listResidue/'
                       '{entry}residue')
            )
        )

        if not all_res:
            print('No mappable residues found', pdb_id, pdb_chain)
            return None

        # make sure to only extract PDBe element and that
        # the PDB and UniProt elements have the right IDs
        pdb_nums = []
        uniprot_nums = []
        for res in all_res:
            if res.attrib['dbSource'] == 'PDBe':
                pdb_record =  None
                uniprot_record = None
                for db_record in list(res):
                    if db_record.attrib['dbSource'] == 'PDB':
                        pdb_record = db_record
                    if db_record.attrib['dbSource'] == 'UniProt':
                        uniprot_record = db_record
                if pdb_record is None or uniprot_record is None:
                    continue
                if pdb_record.attrib['dbResNum'] == 'null' or \
                pdb_record.attrib['dbAccessionId'] != pdb_id or \
                pdb_record.attrib['dbChainId'] != pdb_chain or \
                uniprot_record.attrib['dbAccessionId'] != uniprot_id:
                    continue
                else:
                    pdb_nums.append(pdb_record.attrib['dbResNum'])
                    uniprot_nums.append(uniprot_record.attrib['dbResNum'])

        pdb_to_uniprot = {}
        for x, y in zip(pdb_nums, uniprot_nums):
            try:
                pdb_to_uniprot[int(x)] = int(y)
            except ValueError:
                continue

        return pdb_to_uniprot

    def by_alignment(self, min_overlap=20, reduce_chains=False, **kwargs):
        """
        @TODO

        Parameters
        ----------
        min_overlap
        reduce_chains
        kwargs

        Returns
        -------

        """
        pass


def main():
    """
    This is a test function.

    Returns
    -------

    """
    sifts_uniprot = '/dors/capra_lab/users/lib14/mtr3d/mtr3d/' \
                     'examples/pdb_chain_uniprot.tsv.gz'

    print('SIFTS mapping table file:', sifts_uniprot)

    sifts_table = SIFTS(sifts_uniprot)

    pdb_id = '9rsa'
    pdb_chain = 'A'

    print('PDB chain requested:', pdb_id + pdb_chain)

    uniprot_id = 'P61823'

    pdb_to_uniprot_mapping = sifts_table.pdb_to_uniprot(
        pdb_id, pdb_chain, uniprot_id)

    print('PDB to UniProt mapping:')
    print('=======================')
    print(pdb_to_uniprot_mapping)

    uniprot_to_pdb_mapping = sifts_table.uniprot_to_pdb(
        uniprot_id, pdb_id, pdb_chain)

    print('UniProt to PDB mapping:')
    print('=======================')
    print(uniprot_to_pdb_mapping)


if __name__ == '__main__':
    main()
