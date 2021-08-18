

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
        # get file format 
        self.init_file_format = self._getFileFormat()


        # load OB molecule
        self.ob_mol = None
        self._loadMol()

        # initialize RW mol
        self.rw_mol = None

    def _getFileFormat(self):
        return self.mol_path.split(".")[-1]

    def _loadMol(self):
        # mol = readfile("smi", "myfile.smi").next() # Python 2
        # mol = next(readfile("smi", "myfile.smi"))  # Python 3
        self.ob_mol = next(pybel.readfile(self.init_file_format, self.mol_path))

    def addH(self):
         self.ob_mol.addh()

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
        import sys
        if self.rw_mol is None:
            self.obMol2RWmol()

        file_format = file_path.split(".")[-1]

        if file_format == "xyz":
            rdkit.Chem.rdmolfiles.MolToXYZFile(self.rw_mol, file_path)
        elif file_format == "pdb":
            rdkit.Chem.rdmolfiles.MolToPDBFile(self.rw_mol, file_path)
        else:
            print("Unknown file format")
            sys.exit(1)

    def printTest(self):
        print(len(self._getFileFormat()))

