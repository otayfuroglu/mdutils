

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
        self.mol_path = mol_path

        # load OB molecule
        self.ob_mol = None
        self._loadMol()

        # initialize RW mol
        self.rw_mol = None


        # get file format 
        self.mol_format = self._getFileFormat()


    def _getFileFormat(self):
        return self.mol_path.split(".")[-1]

    def _loadMol(self):
        # mol = readfile("smi", "myfile.smi").next() # Python 2
        # mol = next(readfile("smi", "myfile.smi"))  # Python 3
        self.ob_mol = next(pybel.readfile(self.mol_format, self.mol_path))

    #  def addH(self):
    #       self.ob_mol.addh()

    #  def removeH(self):
    #       self.ob_mol.removeh()

    def _rmFileExist(self, file_path):
        import os
        if os.path.exists(file_path):
            os.remove(file_path)

    def obMol2RWmol(self):
        temp_file_name = "temp_ob_file.mol2"
        self.writeOBMol2File("mol2", temp_file_name)
        self.rw_mol = rdkit.Chem.rdmolfiles.MolFromMol2File(temp_file_name)

        # remove temp
        self._rmFileExist(temp_file_name)

    def writeOBMol2File(self, file_format, file_path):
        with pybel.Outputfile(file_format, file_path) as fl:
            fl.write(self.ob_mol)

    def writeRWMol2File(self, file_path):
        if rw_mol is None:
            self.obMol2RWmol()
        rdkit.Chem.rdmolfiles.MolToPDBFile(rw_mol, file_path)

    def printTest(self):
        print(len(self.ob_mol.atoms))

