

import openbabel
#  from openbabel import pybel
import rdkit


class ligPrep:
    """

    """

    def __init__(self):
        self.mol = None

    def loadMol(self, mol_path):
        self.mol = openbabel.pybel(mol_path)
    
    def test(self):
        pass



