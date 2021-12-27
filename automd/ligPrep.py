

#  import openbabel
#  from openbabel import pybel
import rdkit
from  rdkit import Chem
from  rdkit.Chem import AllChem
import os_util
import sys

class ligPrep:
    """

    """

    def __init__(self, mol_path):
        self.mol_path = mol_path

        # initialize RW mol
        self.rw_mol = None
        self._loadRWMol()

        # initialize calcultor
        self.calculator = None

        # initialize geom opt paprameters
        self.maxiter = None
        self.fmax = None

    def _getFileFormat(self, file_path=None):
        if file_path:
            return file_path.split(".")[-1]

        return self.mol_path.split(".")[-1]

    def _loadRWMol(self):
        if self._getFileFormat() != "mol2":
            print("Error: File fomat is NOT mol2")
            sys.exit(1)
        self.rw_mol = Chem.rdmolfiles.MolFromMol2File(self.mol_path)


    def _addHwithRD(self):
        self.rw_mol = rdkit.Chem.rdmolops.AddHs(self.rw_mol, addCoords=True)

    def _rmFileExist(self, file_path):
        import os
        if os.path.exists(file_path):
            os.remove(file_path)

    def writeRWMol2File(self, file_path):

        # add missing H
        self._addHwithRD()

        file_format = self._getFileFormat(file_path)

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
                         pruneRmsThresh=0.1,
                         mmCalculator=False,
                         optimization=False,
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
            if optimization:
                e = self._geomOptimization(mol, conformerId)
            else:
                if mmCalculator:
                    e = self._calcEnergyWithMM(mol, conformerId, 100)["energy_abs"]
                else:
                    e = self._calcSPEnergy(mol, conformerId)

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


    def setG16Calculator(self, label, chk, xc, basis, scf, multiplicity, extra):
        from ase.calculators.gaussian import Gaussian

        self.calculator = Gaussian(
            label=label,
            chk=chk,
            xc=xc,
            basis=basis,
            scf=scf,
            charge=Chem.rdmolops.GetFormalCharge(self.rw_mol),
            mult=multiplicity,
            extra=extra,
        )

    def setANI2XCalculator(self):
        import torchani
        import torch
        print("Nuber of CUDA devices: ", torch.cuda.device_count())
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

        self.calculator = torchani.models.ANI2x().to(device).ase()

    def _calcSPEnergy(self, mol, conformerId):

        if self.calculator is None:
            print("Error: Calculator not found. Please set any calculator")
            sys.exit(1)

        ase_atoms = self._rwConformer2AseAtoms(mol, conformerId)
        #  from ase.io import write
        #  write("test_ase_atoms.xyz", ase_atoms)
        ase_atoms.set_calculator(self.calculator)

        return ase_atoms.get_potential_energy()

    def calcSPEnergy(self, atoms):

        if self.calculator is None:
            print("Error: Calculator not found. Please set any calculator")
            sys.exit(1)
        atoms.set_calculator(self.calculator)
        return atoms.get_potential_energy()

    def setOptParams(self, fmax, maxiter):
        self.maxiter = maxiter
        self.fmax = fmax

    def _geomOptimization(self, mol, conformerId):
        from ase.optimize import BFGS

        if self.calculator is None:
            print("Error: Calculator not found. Please set any calculator")
            sys.exit(1)


        ase_atoms = self._rwConformer2AseAtoms(mol, conformerId)
        #  from ase.io import write
        #  write("test_ase_atoms.xyz", ase_atoms)

        if self.fmax is None or self.maxiter is None:
            print("Error setting geometry optimizatian parameters for ASE. Please do it")
            exit(1)


        ase_atoms.set_calculator(self.calculator)
        dyn = BFGS(ase_atoms)
        dyn.run(fmax=self.fmax,steps=self.maxiter)

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
