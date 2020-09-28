import numpy as np

RESIDUE_POLARITY_TYPE = {
    'PHE': 0,
    'MET': 0,
    'ILE': 0,
    'LEU': 0,
    'VAL': 0,
    'CYS': 0,
    'PRO': 0,
    'ALA': 0,
    'GLY': 0,
    'THR': 0.5,
    'TRP': 0.5,
    'SER': 0.5,
    'TYR': 0.5,
    'GLN': 0.5,
    'ASN': 0.5,
    'HIS': 1,
    'GLU': -1,
    'LYS': 1,
    'ASP': -1,
    'ARG': 1
}

BACKBONE_ATOMS = {'N', 'CA', 'C', 'O'}


class Contact:
    """
    """

    def __init__(self, res_a, res_b):
        """
        Parameters
        ----------
        res_a : Residue
            The first residue making this contact. A Biopython Residue object.
        res_b : Residue
            The second residue making this contact. A Biopython Residue object.
        """
        self._res_a = res_a
        self._res_b = res_b

    def __str__(self):
        """
        Returns
        -------
        str
            A user friendly representation of this contact.
        """

    def get_res_a(self):
        """
        Returns
        -------
        Residue
            The first residue of this contact.
        """
        return self._res_a

    def get_res_b(self):
        """
        Returns
        -------
        Residue
            The second residue of this contact.
        """
        return self._res_b

    def get_res_a_index(self):
        """
        Returns
        -------
        int
            Sequence index of the first residue.
        -------

        """
        return self._res_a.get_id()[1]

    def get_res_b_index(self):
        """
        Returns
        -------
        int
            Sequence index of the second residue.
        """
        return self._res_b.get_id()[1]

    def get_shortest_distance(self):
        """
        Returns
        -------
        float
            The smallest of the distances between the heavy atoms
            of residue a and those of residue b.
        """
        distances = [atom_a - atom_b for atom_a in self._res_a.get_list()
                     for atom_b in self._res_b.get_list()]
        return min(distances)

    def get_cb_distance(self):
        """
        Returns
        -------
        float
            The distance between the Cb atoms of the contacting residues.
        """
        try:
            res_a_cb = self._res_a['CB']
        except KeyError:
            res_a_cb = self._res_a['CA']
        try:
            res_b_cb = self._res_b['CB']
        except KeyError:
            res_b_cb = self._res_b['CA']
        return res_a_cb - res_b_cb

    @staticmethod
    def get_sidechain_centroid(res):
        """
            C = (v1 + ... + vn) / n
        """
        # return CA coordinates if residue is GLY
        if res.get_resname() == 'GLY':
            return res['CA'].get_coord()
        # for other residues, return sidechain centroid
        centroid = np.array([0.0, 0.0, 0.0])
        number_sidechain_atoms = 0
        for atom in res.get_list():
            if atom.get_name() not in BACKBONE_ATOMS:
                centroid += atom.get_coord()
                number_sidechain_atoms += 1
        return centroid / number_sidechain_atoms

    def get_centroid_distance(self):
        """
        Returns
        -------
        float
            The distance between the side-chain centroids of the contacting
            residues.
        """
        centroid_a = self.get_sidechain_centroid(self._res_a)
        centroid_b = self.get_sidechain_centroid(self._res_b)
        return np.linal.norm(centroid_a - centroid_b)

    def is_electrostatic(self):
        """
        Determines whether the contact is of is_electrostatic nature.

        Returns
        -------
        bool
            True if the contacting residues are oppositely charged.
        """
        res_a_type = RESIDUE_POLARITY_TYPE[self._res_a.get_resname()]
        res_b_type = RESIDUE_POLARITY_TYPE[self._res_b.get_resname()]
        return res_a_type * res_b_type < 0

    def is_nonpolar(self):
        """
        Determines whether the contact is of nonpolar nature.

        Returns
        -------
        bool
            True if none of the contact residues is polar.
        """
        res_a_type = RESIDUE_POLARITY_TYPE[self._res_a.get_resname()]
        res_b_type = RESIDUE_POLARITY_TYPE[self._res_b.get_resname()]
        return res_a_type == 0 and res_b_type == 0

    def is_polar(self):
        """
        Determines whether the contact is of polar nature.

        Returns
        -------
        bool
            True if the either of the contact residues is polar.
        """
        return not self.is_nonpolar()

    def is_local(self):
        """
        @TODO
        """
        pass

    def is_global(self):
        """
        Determines whether a contact is global.

        Returns
        -------
        bool
            True if the contact is global.
        """
        return not self.is_local()

