#
from openbabel import openbabel
from ase import Atoms
from ase.io import write

import rdkit
from  rdkit import Chem
from  rdkit.Chem import AllChem
#  import os_util
from collections import defaultdict
import sys, os
import numpy as np
import pandas as pd

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
        try:
            self._loadMolWithRW(tmp_file_name)
        except:
            self._loadMolWithRW(tmp_file_name, sanitize=False)

            # correction  kekuleize error (especially N in aromatic ring)
            print("\nWarning!: There is kekulize error, ingnored sanitize and kekulized for N atom which is in aromatic ring\n")
            for i, atom in enumerate(self.rw_mol.GetAtoms()):
                if atom.GetSymbol() == "N" and atom.GetIsAromatic():
                    print("Aromatic N atom idex: ",i+1)
                    self.rw_mol.GetAtomWithIdx(i+1).SetNumExplicitHs(1)
        # remove temp
        self._rmFileExist(tmp_file_name)

        # NOTE
        #  ase_atoms = self.obMol2AseAtoms(ob_mol)

        # optmization for just added H
        if self.addH:
            self.setOptMethod(opt_method="LBFGS")
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

    def _getNumConfs(self, scaled=1):
        n_torsions = len(self._getTorsionPoints())
        if n_torsions > 10:
            return 5000
        else:
            return int(scaled * 2 ** n_torsions)

    def _getConfDistMatrix(self, mol, conformerIds):

        #  conf_list = []
        #  for conformerId in conformerIds:
        #      print(conformerId)
        #      conf = mol.GetConformer(conformerId)
        #      #  positions = conf.GetPositions()
        #      #  print(positions)
        #      conf.SetProp("_Name", str(conformerId))
        #      conf_list.append(conf)

        #  print(conf_list)
        n_mol=len(conformerIds)
        if n_mol <= 1:
            print("Clustering do not applied.. There is just one conformer")
            return 0
        conf_dist_matrix = np.empty(shape=(n_mol, n_mol))
        for conf1_id in conformerIds:
            for conf2_id in conformerIds:
                conf_dist_matrix[conf1_id, conf2_id] = AllChem.GetConformerRMS(mol, conf1_id, conf2_id, prealigned=False)

        return conf_dist_matrix

    def _getClusterKmeansFromConfIds(self, conformerIds, dist_matrix, n_group):
        from scipy.cluster.vq import kmeans, vq, whiten

        cluster_conf_id = defaultdict(list)
        whitened = whiten(dist_matrix)
        centroids, _ = kmeans(whitened, n_group)
        cluster, _ = vq(whitened,centroids)
        for key, value in zip(cluster, conformerIds):
            cluster_conf_id[key].append(value)

        return cluster_conf_id

    def _getClusterRMSDFromFiles(self, conf_dir, rmsd_thresh):
        from scipy.cluster.hierarchy import linkage, fcluster
        suppl_list = {Chem.SDMolSupplier(f"{conf_dir}/{fl_name}"):
                      fl_name for fl_name in os.listdir(conf_dir)
                      if fl_name.endswith(".sdf")}
        mol_list = []
        for suppl, fl_name in suppl_list.items():
            mol = next(suppl)
            mol.SetProp("_Name", fl_name)
            mol_list.append(mol)

        n_mol=len(mol_list)
        if n_mol <= 1:
            print("Clustering do not applied.. There is just one conformer")
            return 0
        dist_matrix=np.empty(shape=(n_mol, n_mol))
        for i, mol1 in enumerate(mol_list):
            for j, mol2 in enumerate(mol_list):
                # first aling each conformer and then calculate rmsd
                #  Chem.rdMolAlign.AlignMol(mol1, mol2)
                dist_matrix[i, j] = Chem.rdMolAlign.GetBestRMS(mol1, mol2)
        linked = linkage(dist_matrix,'complete')
        cluster_conf = defaultdict(list)
        labelList = [mol.GetProp('_Name') for mol in mol_list]
        for key, fl_name in zip(fcluster(linked, rmsd_thresh, criterion='distance'), labelList):
            cluster_conf[key].append(fl_name)

            # save clusturedd files seperately
            directory = f"{conf_dir}/cluster_{key}"
            if not os.path.exists(directory):
                os.mkdir(directory)
            file_path = f"{directory}/{fl_name}"
            for mol in mol_list:
                if mol.GetProp('_Name') == fl_name:
                    mol = mol
                    break
            with Chem.rdmolfiles.SDWriter(file_path) as writer:
                writer.write(mol)

        return cluster_conf

    def _pruneOptConfs(self, cluster_conf, confs_energies, conf_dir):
        i = 0
        for fl_names in cluster_conf.values():
            for j, fl_name in enumerate(fl_names):
                e = float(confs_energies.loc[confs_energies["FileName"] == fl_name, " Energy(eV)"].item())
                if i == 0:
                    global_minE = e
                    global_minE_file = fl_name
                else:
                    if global_minE > e:
                        global_minE = e
                        global_minE_file = fl_name
                i += 1

                if j == 0:
                    minE = e
                    minE_file = fl_name
                else:
                    if minE > e:
                        minE = e
                        minE_file = fl_name
            fl_names.remove(minE_file)
            os.rename(f"{conf_dir}/{minE_file}", f"{conf_dir}/pruned_{minE_file}")
            if len (fl_names) != 0:
                for rm_file in fl_names:
                    print("Removed", rm_file)
                    os.remove(f"{conf_dir}/{rm_file}")

        # file which has global minimum enery renamed 
        os.rename(f"{conf_dir}/pruned_{global_minE_file}", f"{conf_dir}/global_minE_{global_minE_file}")

    def genMinEGonformer(self, file_path,
                         numConfs=100,
                         ETKDG=False,
                         maxAttempts=10000,
                         pruneRmsThresh=0.1,
                         mmCalculator=False,
                         optimization_conf=False,
                         opt_prune_rms_thresh=1.0,
                         saveConfs=True,
                         useExpTorsionAnglePrefs=True,
                         useBasicKnowledge=True,
                         enforceChirality=True
                        ):

        import copy

        #  self.addHwithRD()
        mol = copy.deepcopy(self.rw_mol)
        if numConfs is 0 or numConfs < self._getNumConfs(scaled=10):
            numConfs = self._getNumConfs(scaled=10)
            print(f"Maximum number of conformers setting to {numConfs}")

        if ETKDG:
            confs = rdkit.Chem.AllChem.EmbedMultipleConfs(
                mol, numConfs=numConfs)
        else:
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

        # file for saving energies
        file_csv = open("%s/all_confs_sp_energies.csv" %self.WORK_DIR, "w")
        print("FileName, Energy(eV)", file=file_csv)

        print("Number of generated conformation: %d" %len(confs))

        dist_matrix = self._getConfDistMatrix(mol, confs)

        cluster_conf_id = self._getClusterKmeansFromConfIds(confs, dist_matrix,
                                           n_group=self._getNumConfs(scaled=1)
                                          )

        #  for k-means clutering
        minEConformerIDs = []
        for cluster, clustered_confIds in cluster_conf_id.items():

            if saveConfs:
                CONF_DIR = self.WORK_DIR + f"/confs_cluster_{cluster}"
                if not os.path.exists(CONF_DIR):
                    os.mkdir(CONF_DIR)

            for i, conformerId  in enumerate(clustered_confIds):
                #  print("%d. conformer processing..." %i)
                if saveConfs:
                    prefix = ""
                    conf_file_path = "%s/conf_%d.sdf"%(CONF_DIR, conformerId)
                    self._writeConf2File(mol, conformerId, conf_file_path)

                #create ase atoms
                ase_atoms = self._rwConformer2AseAtoms(mol, conformerId)
                if mmCalculator:
                    e = self._calcEnergyWithMM(mol, conformerId, 100)["energy_abs"]
                else:
                    e = self._calcSPEnergy(mol, conformerId)

                if i == 0:
                    minE = e
                    minEConformerID = conformerId
                    minE_ase_atoms = ase_atoms
                else:
                    if minE > e:
                        minE = e
                        minEConformerID = conformerId
                        minE_ase_atoms = ase_atoms
                print("%sconf_%d.sdf, %s"%(prefix, conformerId, e), file=file_csv)

            minEConformerIDs.append(minEConformerID)

            #  print(minE, minEGonformerID)
            #  if numConfs > 1:
            #  if not optimization_conf:
            #      _, ase_atoms = self._geomOptimizationConf(mol, minEGonformerID)
            #  else:
            #      ase_atoms = self._rwConformer2AseAtoms(mol, minEGonformerID)

            # compeare initial struct and minE conf. energies
            #  if self.calcSPEnergy() >  minE:
            #      # assign minE conformer to self rw mol
            #      self.rw_mol = self.aseAtoms2rwMol(minE_ase_atoms)
            #      print("Selected minimun energy conformer")
            #  else:
            #      print("!!! INITIAL structure as is minimun energy conformer !!!")

        # test
        assert len(minEConformerIDs) == len(cluster_conf_id.keys())
        # close to csv file
        file_csv.close()

        MIN_E_CONF_DIR = self.WORK_DIR + "/opt_minE_confs"
        if not os.path.exists(MIN_E_CONF_DIR):
            os.mkdir(MIN_E_CONF_DIR)

        if optimization_conf:
            opt_file_csv = open("%s/opt_confs_energies.csv" %MIN_E_CONF_DIR, "w")
            print("FileName, Energy(eV)", file=opt_file_csv)

            for i, conformerId  in enumerate(minEConformerIDs):
                e, ase_atoms = self._geomOptimizationConf(mol, conformerId)
                prefix = "opt_"
                conf_file_path = "%s/%sconf_%d.xyz"%(MIN_E_CONF_DIR, prefix, conformerId)

                # save optimized structure  with ase
                write(conf_file_path, ase_atoms)

                # save optimized structure  with rdkit as sdf
                #  with Chem.rdmolfiles.SDWriter(conf_file_path) as writer:
                #      writer.write(self.aseAtoms2rwMol(ase_atoms))

                print("%sconf_%d.sdf, %s"%(prefix, conformerId, e), file=opt_file_csv)
            opt_file_csv.close()

            # cluster and prune opitimzed confs by RMSD
            confs_energies = pd.read_csv(f"{MIN_E_CONF_DIR}/opt_confs_energies.csv")
            print(confs_energies)
            cluster_conf = self._getClusterRMSDFromFiles(MIN_E_CONF_DIR, rmsd_thresh=opt_prune_rms_thresh)
            if cluster_conf != 0:
                self._pruneOptConfs(cluster_conf, confs_energies, MIN_E_CONF_DIR)

                #  #create ase atoms
                #  ase_atoms = self._rwConformer2AseAtoms(mol, conformerId)
                #  if mmCalculator:
                #      e = self._calcEnergyWithMM(mol, conformerId, 100)["energy_abs"]
                #  else:
                #      e = self._calcSPEnergy(mol, conformerId)

                #  if i == 0:
                #      minE = e
                #      minEConformerID = conformerId
                #      minE_ase_atoms = ase_atoms
                #  else:
                #      if minE > e:
                #          minE = e
                #          minEConformerID = conformerId
                #          minE_ase_atoms = ase_atoms


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
        if self.opt_method is None or self.opt_method=="lbfgs":
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

        file_csv = open("%s/optmized_energy.csv" %self.WORK_DIR, "w")
        print(ase_atoms.get_potential_energy(), ",eV", file=file_csv)

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
