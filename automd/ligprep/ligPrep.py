#
from openbabel import openbabel
from ase import Atoms
from ase.io import write

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

        # for activete g16 optmization algorithm
        self.optG16 = False

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

        # initialize optimization method
        self.opt_method = None


    def _getFileFormat(self, file_path=None):
        if file_path:
            return file_path.split(".")[-1]

        return self.mol_path.split(".")[-1]

    def _loadRWMol(self):
        #  if self._getFileFormat() == "mol2" and not self.addH:
        #      self._loadMolWithRW(self.mol_path)
        #  else:
        self._loadMolWithOB()

    def _loadMolWithRW(self, mol_path, sanitize=True):
        rd_mol = Chem.rdmolfiles.MolFromMol2File(mol_path, sanitize=sanitize, removeHs=False)
        self.rw_mol = Chem.RWMol(rd_mol)

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

        #  # openbabel file to rdkit mol2 file
        tmp_file_name = "tmp_ob_file.mol2"
        obConversion.WriteFile(ob_mol, tmp_file_name)

        # laod as RW file
        self._loadMolWithRW(tmp_file_name)
        # remove temp
        self._rmFileExist(tmp_file_name)

        # NOTE
        #  ase_atoms = self.obMol2AseAtoms(ob_mol)

        # optmization for just added H
        if self.addH:
            self.setOptParams(fmax=0.05, maxiter=200)
            self.setANI2XCalculator()
            self.geomOptimization(fix_heavy_atoms=True)

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

    def _getTorsionPoints(self):
        from rdkit.Chem import TorsionFingerprints
        torsion_points = []
        for torsions_list in TorsionFingerprints.CalculateTorsionLists(self.rw_mol):
            for torsions in torsions_list:
                if 180 in torsions:
                    torsion_points.append(torsions[0][0])
        return(torsion_points)

    def _getNunConfs(self):
        n_torsions = len(self._getTorsionPoints())
        if n_torsions > 10:
            return 10
        else:
            return 2 ** n_torsions


    def genMinEGonformer(self, file_path,
                         numConfs=0,
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
        if numConfs is 0:
            numConfs = self._getNunConfs()
            print(f"Maximum number of conformers setting to {numConfs}")


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

        CONF_DIR = self.WORK_DIR + "/conformers"
        if not os.path.exists(CONF_DIR):
            os.mkdir(CONF_DIR)

        # file for saving energies
        file_csv = open("%s/confs_energies.csv" %CONF_DIR, "w")
        print("FileName, Energy(eV)", file=file_csv)

        print("Number of generated conformation: %d" %len(confs))
        for i, conformerId  in enumerate(confs):
            print("%d. conformer processing..." %i)
            if saveConfs:
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
            print("%sconf_%d.sdf, %s"%(prefix, i, e), file=file_csv)
        file_csv.close()
        #  print(minE, minEGonformerID)
        #  if numConfs > 1:
        #  if not optimization_conf:
        #      _, ase_atoms = self._geomOptimizationConf(mol, minEGonformerID)
        #  else:
        #      ase_atoms = self._rwConformer2AseAtoms(mol, minEGonformerID)

        # compeare initial struct and minE conf. energies
        if self.calcSPEnergy() >  minE:
            # assign minE conformer to self rw mol
            self.rw_mol = self.aseAtoms2rwMol(minE_ase_atoms)
            print("Selected minimun energy conformer")
        else:
            print("!!! INITIAL structure as is minimun energy conformer !!!")

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

    def setG16Calculator(self, label, chk, nprocs, xc, basis, scf, addsec=None, extra=None):
        from ase.calculators.gaussian import Gaussian
        self.optG16 = True

        self.calculator = Gaussian(
            label=label,
            chk=chk,
            nprocshared=nprocs,
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

    def calcSPEnergy(self):

        if self.calculator is None:
            print("Error: Calculator not found. Please set any calculator")
            sys.exit(1)
        ase_atoms= self.rwMol2AseAtoms()
        ase_atoms.set_calculator(self.calculator)
        return ase_atoms.get_potential_energy()

    def setOptParams(self, fmax, maxiter):
        self.maxiter = maxiter
        self.fmax = fmax

    def setOptMethod(self, opt_method):
        self.opt_method = opt_method.lower()

    def _getOptMethod(self, ase_atoms):
        if self.opt_method is None or self.opt_method=="lfbgs":
            from ase.optimize import LBFGS
            return LBFGS(ase_atoms)
        elif self.opt_method=="bfgs":
            from ase.optimize import BFGS
            return BFGS(ase_atoms)
        elif self.opt_method=="fire":
            from ase.optimize import FIRE
            return FIRE(ase_atoms)
        elif self.opt_method=="gpmin":
            from ase.optimize import GPMin
            return GPMin(ase_atoms)
        elif self.opt_method=="berny":
            from ase.optimize import Berny
            return Berny(ase_atoms)

    def _geomOptimizationConf(self, mol, conformerId):
        from ase.calculators.gaussian import GaussianOptimizer, Gaussian

        if self.calculator is None:
            print("Error: Calculator not found. Please set any calculator")
            sys.exit(1)


        ase_atoms = self._rwConformer2AseAtoms(mol, conformerId)
        #  from ase.io import write
        #  write("test_ase_atoms.xyz", ase_atoms)

        if self.fmax is None or self.maxiter is None:
            print("Error setting geometry optimizatian parameters for ASE. Please do it")
            exit(1)

        if self.optG16:
            dyn =  GaussianOptimizer(ase_atoms, self.calculator)
            dyn.run(fmax='tight', steps=self.maxiter)
        else:
            ase_atoms.set_calculator(self.calculator)
            dyn = self._getOptMethod(ase_atoms)
            dyn.run(fmax=self.fmax, steps=self.maxiter)

        #  ase_atoms.set_calculator(self.calculator)
        #  dyn = LBFGS(ase_atoms)
        #  dyn.run(fmax=self.fmax,steps=self.maxiter)

        return ase_atoms.get_potential_energy(), ase_atoms

    def geomOptimization(self, fix_heavy_atoms=False):
        from ase.calculators.gaussian import GaussianOptimizer

        if self.calculator is None:
            print("Error: Calculator not found. Please set any calculator")
            sys.exit(1)
        if self.fmax is None or self.maxiter is None:
            print("Error setting geometry optimizatian parameters for ASE. Please do it")
            exit(1)

        ase_atoms = self.rwMol2AseAtoms()
        if fix_heavy_atoms:
            from ase.constraints import FixAtoms
            c = FixAtoms(indices=[atom.index for atom in ase_atoms if atom.symbol != 'H'])
            ase_atoms.set_constraint(c)

        if self.optG16:
            dyn =  GaussianOptimizer(ase_atoms, self.calculator)
            dyn.run(fmax='tight', steps=self.maxiter)
        else:
            ase_atoms.set_calculator(self.calculator)
            #  self.dyn = LBFGS(ase_atoms)
            dyn = self._getOptMethod(ase_atoms)
            dyn.run(fmax=self.fmax, steps=self.maxiter)

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

    def obMol2AseAtoms(self, ob_mol):
        from ase import Atom
        ase_atoms = Atoms()
        for i in range(ob_mol.NumAtoms()):
            obatom = ob_mol.GetAtom(i + 1)
            ase_atoms.append(Atom(obatom.GetAtomicNum(),
                              [obatom.GetX(),
                               obatom.GetY(),
                               obatom.GetZ()]
                             ))
        return ase_atoms

    def aseAtoms2rwMol(self, ase_atoms):

        #  from aseAtoms2rdMol import fromASE, to_rdmol, ase2rdmol
        write("tmp.pdb", ase_atoms)
        rd_mol = Chem.rdmolfiles.MolFromPDBFile("tmp.pdb", removeHs=False)
        self._rmFileExist("tmp.pdb")

        return AllChem.AssignBondOrdersFromTemplate(self.rw_mol, rd_mol)


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
