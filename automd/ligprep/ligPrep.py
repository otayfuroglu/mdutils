#
from openbabel import openbabel
from ase import Atoms
from ase.io import write
from ase.optimize import BFGS

import rdkit
from  rdkit import Chem
from  rdkit.Chem import AllChem
#  import os_util
import sys, os

class ligPrep:
    """

    """

    def __init__(self, mol_path, addH, WORK_DIR):
        self.mol_path = mol_path
        self.WORK_DIR = WORK_DIR

        # set add missing H
        self.addH = addH

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
        if self._getFileFormat() == "mol2" and not self.addH:
            self._loadMolWithRW(self.mol_path)
        else:
            self._loadMolWithOB()


    def _loadMolWithRW(self, mol_path):
        self.rw_mol = Chem.rdmolfiles.MolFromMol2File(mol_path, removeHs=False)

    def _loadMolWithOB(self):

        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats(self._getFileFormat(), "mol2")
        ob_mol = openbabel.OBMol()
        obConversion.ReadFile(ob_mol, self.mol_path)

        #add hydrogen with openbabel
        if self.addH:
            if self._getFileFormat() == "xyz":
                print("Error: Cannot add hydrogen atoms to XYZ file format!!!")
                exit(1)
            ob_mol.AddHydrogens()

        # openbabel file to rdkit mol2 file
        tmp_file_name = "tmp_ob_file.mol2"
        obConversion.WriteFile(ob_mol, tmp_file_name)

        # laod as RW file
        self._loadMolWithRW(tmp_file_name)

        # remove temp
        self._rmFileExist(tmp_file_name)

    def addHwithRD(self):
        self.rw_mol = rdkit.Chem.rdmolops.AddHs(self.rw_mol, addCoords=True)

    def _rmFileExist(self, file_path):
        if os.path.exists(file_path):
            os.remove(file_path)

    def writeRWMol2File(self, file_path):

        # add missing H
        #  self.addHwithRD()

        file_format = self._getFileFormat(file_path)

        if file_format == "xyz":
            rdkit.Chem.rdmolfiles.MolToXYZFile(self.rw_mol, file_path)
        elif file_format == "pdb":
            rdkit.Chem.rdmolfiles.MolToPDBFile(self.rw_mol, file_path)
        elif file_format == "sdf":
            with Chem.rdmolfiles.SDWriter(file_path) as writer:
                writer.write(self.rw_mol)

        else:
            print("Unknown file format")
            sys.exit(1)

    def _writeConf2File(self, mol, conformerId, file_path):
        with rdkit.Chem.SDWriter(file_path) as w:
            w.write(mol, conformerId)
            w.flush()
            w.close()

    def genMinEGonformer(self, file_path,
                         numConfs=20,
                         maxAttempts=1000,
                         pruneRmsThresh=0.1,
                         mmCalculator=False,
                         optimization_conf=False,
                         saveConfs=True,
                         useExpTorsionAnglePrefs=True,
                         useBasicKnowledge=True,
                         enforceChirality=True):

        import copy

        #  self.addHwithRD()
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

        print("Number of generated conformation: %d" %len(confs))
        for i, conformerId  in enumerate(confs):
            print("%d. conformer processing..." %i)
            if saveConfs:
                CONF_DIR = self.WORK_DIR + "/conformers"
                if not os.path.exists(CONF_DIR):
                    os.mkdir(CONF_DIR)

                prefix = ""
                if optimization_conf:
                    prefix = "opt_"
                conf_file_path = "%s/%sconf_%d.sdf"%(CONF_DIR, prefix, i)
                self._writeConf2File(mol, conformerId, conf_file_path)

            if optimization_conf:
                e, ase_atoms = self._geomOptimizationConf(mol, conformerId)
            else:
                #create ase atoms
                ase_atoms = self._rwConformer2AseAtoms(mol, conformerId)
                if mmCalculator:
                    e = self._calcEnergyWithMM(mol, conformerId, 100)["energy_abs"]
                else:
                    e = self._calcSPEnergy(mol, conformerId)

            if i == 0:
                minE = e
                #  minEGonformerID = conformerId
                minE_ase_atoms = ase_atoms
            else:
                if minE > e:
                    minE = e
                    #  minEGonformerID = conformerId
                    minE_ase_atoms = ase_atoms
        #  print(minE, minEGonformerID)
        #  if numConfs > 1:
        #  if not optimization_conf:
        #      _, ase_atoms = self._geomOptimizationConf(mol, minEGonformerID)
        #  else:
        #      ase_atoms = self._rwConformer2AseAtoms(mol, minEGonformerID)

        # assign minE conformer to self rw mol
        self._rw_mol = self.aseAtoms2rwMol(minE_ase_atoms)

        return minE_ase_atoms

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

    def setG16Calculator(self, label, chk, xc, basis, scf, addsec, extra):
        from ase.calculators.gaussian import Gaussian

        self.calculator = Gaussian(
            label=label,
            chk=chk,
            xc=xc,
            basis=basis,
            scf=scf,
            addsec=addsec,
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

    def _geomOptimizationConf(self, mol, conformerId):

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

        return ase_atoms.get_potential_energy(), ase_atoms

    def geomOptimization(self, ase_atoms):

        if self.calculator is None:
            print("Error: Calculator not found. Please set any calculator")
            sys.exit(1)
        if self.fmax is None or self.maxiter is None:
            print("Error setting geometry optimizatian parameters for ASE. Please do it")
            exit(1)

        ase_atoms.set_calculator(self.calculator)
        dyn = BFGS(ase_atoms)
        dyn.run(fmax=self.fmax,steps=self.maxiter)

        self.rw_mol = self.aseAtoms2rwMol(ase_atoms)


        #  return ase_atoms.get_potential_energy(), ase_atoms

    def _rwConformer2AseAtoms(self, mol, conformerId):

        mol = mol.GetConformer(conformerId)

        atom_species = [atom.GetAtomicNum() for atom in mol.GetOwningMol().GetAtoms()]
        positions = mol.GetPositions()

        return Atoms(atom_species, positions)

    def rwMol2AseAtoms(self):

        atom_species = [atom.GetAtomicNum() for atom in self.rw_mol.GetAtoms()]

        conf = self.rw_mol.GetConformer()
        positions = [conf.GetAtomPosition(i) for i in range(len(atom_species))]
        #  positions = self.rw_mol.GetConformers()[0].GetPositions()

        return Atoms(atom_species, positions)

    def aseAtoms2rwMol(self, ase_atoms):

        write("tmp.pdb", ase_atoms)
        rw_mol = Chem.rdmolfiles.MolFromPDBFile("tmp.pdb", removeHs=False)
        self._rmFileExist("tmp.pdb")
        return rw_mol

    def writeAseAtoms(self, file_path):
        ase_atoms = self.rwMol2AseAtoms()

        # write mol to xyz file by ase
        write(file_path, ase_atoms)


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
