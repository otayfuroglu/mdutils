

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
        self._obMol2RWmol()

    def _getFileFormat(self):
        return self.mol_path.split(".")[-1]

    def _loadMol(self):
        # mol = readfile("smi", "myfile.smi").next() # Python 2
        # mol = next(readfile("smi", "myfile.smi"))  # Python 3
        self.ob_mol = next(pybel.readfile(self.init_file_format, self.mol_path))

    def addHwithOB(self):
         self.ob_mol.addh()

    def _addHwithRD(self):
        self.rw_mol = rdkit.Chem.rdmolops.AddHs(self.rw_mol, addCoords=True)

    def _rmFileExist(self, file_path):
        import os
        if os.path.exists(file_path):
            os.remove(file_path)

    def _obMol2RWmol(self):
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

        # add missing H
        self._addHwithRD()

        file_format = file_path.split(".")[-1]

        if file_format == "xyz":
            rdkit.Chem.rdmolfiles.MolToXYZFile(self.rw_mol, file_path)
        elif file_format == "pdb":
            rdkit.Chem.rdmolfiles.MolToPDBFile(self.rw_mol, file_path)
        else:
            print("Unknown file format")
            sys.exit(1)

    def genMinEGonformer(self, file_path,
                         numConfs=20,
                         maxAttempts=1000,
                         pruneRmsThresh=0.5,
                         useExpTorsionAnglePrefs=True,
                         useBasicKnowledge=True,
                         enforceChirality=True):

        import copy

        self._addHwithRD()
        mol = copy.deepcopy(self.rw_mol)

        confs = rdkit.Chem.AllChem.EmbedMultipleConfs(
            mol,
            numConfs=numConfs,
            maxAttempts=maxAttempts,
            pruneRmsThresh=pruneRmsThresh,
            useExpTorsionAnglePrefs=useExpTorsionAnglePrefs,
            useBasicKnowledge=useBasicKnowledge,
            enforceChirality=enforceChirality,
            numThreads=0,
        )

        for i, conformerId  in enumerate(confs):
            #e = self._calcEnergyWithMM(mol, conformerId, 100)["energy_abs"]
            e = self._gemOptWithG16(mol, conformerId, 'B3LYP', "sto-3g", fmax=0.05)
            #e = self._calcEnergyWithG16(mol, conformerId, 'B3LYP', "sto-3g")
            if i == 0:
                minE = e
                minEGonformerID = conformerId
            else:
                if minE > e:
                    minE = e
                    minEGonformerID = conformerId
        #  print(minE, minEGonformerID)

        with rdkit.Chem.SDWriter(file_path) as w:
            w.write(mol, minEGonformerID)
            w.flush()
            w.close()

    def _calcEnergyWithMM(self, mol, conformerId, minimizeIts):
        ff = rdkit.Chem.AllChem.MMFFGetMoleculeForceField(
            mol,
            rdkit.Chem.AllChem.MMFFGetMoleculeProperties(mol),
            confId=conformerId)
        ff.Initialize()
        ff.CalcEnergy()
        results = {}
        if minimizeIts > 0:
            results["converged"] = ff.Minimize(maxIts=minimizeIts)
        results["energy_abs"] = ff.CalcEnergy()
        return results

    def _calcEnergyWithG16(self, mol, conformerId, xc, basiset):
        from ase.calculators.gaussian import Gaussian

        ase_atoms = self._rwConformer2AseAtoms(mol, conformerId)
        #  from ase.io import write
        #  write("test_ase_atoms.xyz", ase_atoms)

        ase_atoms.set_calculator(
            Gaussian(
                label = "tempG16/gaussionSPE",
                xc=xc,
                basis=basiset,
                scf="maxcycle=100",
                #cpu="0-15",
                nprocshared="32"

            )
        )

        return ase_atoms.get_potential_energy()

    def _gemOptWithG16(self, mol, conformerId, xc, basiset, fmax=0.05):
        from ase.calculators.gaussian import Gaussian
        from ase.optimize import BFGS

        ase_atoms = self._rwConformer2AseAtoms(mol, conformerId)
        #  from ase.io import write
        #  write("test_ase_atoms.xyz", ase_atoms)

        ase_atoms.set_calculator(
            Gaussian(
                label = "tempG16/gaussionSPE",
                xc=xc,
                basis=basiset,
                scf="maxcycle=100",
                #cpu="0-15",
                nprocshared="32"

            )
        )

        dyn = BFGS(ase_atoms)
        dyn.run(fmax)

        return ase_atoms.get_potential_energy()

    def _rwConformer2AseAtoms(self, mol, conformerId):
        from ase import Atoms

        mol = mol.GetConformer(conformerId)

        atom_species = [atom.GetAtomicNum() for atom in mol.GetOwningMol().GetAtoms()]
        positions = mol.GetPositions()

        return Atoms(atom_species, positions)

#      def alignMols(self, mol):
#          rdkit.Chem.rdMolAlign.AlignMol(mol, self.rw_mol)


#  def alignConfs(mol, clust_ids):
#      rmslist = []
#      AllChem.AlignMolConformers(mol, confIds=clust_ids, RMSlist=rmslist)
#      return rmslist
#
#  def calcRMS(mol, ref_mol):
#      rms = Chem.rdMolAlign.CalcRMS(mol, ref_mol)
#      return rms

#  def writeConf2sdf(mol, filename, confId):
#      w = Chem.SDWriter(filename)
#      w.write(mol, confId)
#      w.close()
#
#  def getClusterConf(mol, mode="RMSD", threshold=2.0):
#      if mode == "TFD":
#          dmat = TorsionFingerprints.GetTFDMatrix(mol)
#      else:
#          dmat = AllChem.GetConformerRMSMatrix(mol, prealigned=False)
#      rms_clusters = Butina.ClusterData(dmat, mol.GetNumConformers(), threshold, isDistData=True, reordering=True)
#      return rms_clusters
#
#  #  mol = Chem.MolFromMolFile("./mutemel/article_No1_fin.mol")
#  #  ref_mol = Chem.MolFromMolFile("./mutemel/article_No1_fin_conf_20.sdf")
#  #  print(calcRMS(mol, ref_mol))
