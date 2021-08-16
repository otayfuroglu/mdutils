

import openbabel
from openbabel import pybel
import rdkit
import rdkit.Chem
import rdkit.Chem.AllChem
import os_util

class ligPrep:
    """

    """

    def __init__(self, mol_path):
        self.mol = None
        self.mol_path = mol_path


        # get file format 
        self.mol_format = self._getFileFormat()

        # load molecule
        self._loadMol()

    def _getFileFormat(self):
        return self.mol_path.split(".")[-1]

    def _loadMol(self):
        # mol = readfile("smi", "myfile.smi").next() # Python 2
        # mol = next(readfile("smi", "myfile.smi"))  # Python 3
        self.mol = next(pybel.readfile(self.mol_format, self.mol_path))

    def addH(self):
         self.mol.addh()

    def removeH(self):
         self.mol.removeh()

    def writeOBMol2File(self, file_format, file_path):
        with pybel.Outputfile(file_format, file_path) as fl:
            fl.write(self.mol)


    def writeRWMol2File(self, file_path):
       rwmol = self.getObmolFromRWmol()
       rdkit.Chem.rdmolfiles.MolToPDBFile(rwmol, file_path)

    def printTest(self):
        print(len(self.mol.atoms))

