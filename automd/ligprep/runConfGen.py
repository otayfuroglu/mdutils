
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import numpy as np
from collections import defaultdict
import pandas as pd


from ase.io import read
from ligPrep import ligPrep
import argparse
import os, sys, shutil
import multiprocessing

nprocs_all = int(multiprocessing.cpu_count())



parser = argparse.ArgumentParser(description="Give something ...")
parser.add_argument("structure_dir", type=str)
parser.add_argument("add_hydrogen", nargs="?", default="No") # args for bool
parser.add_argument("calculator_type", type=str)
parser.add_argument("optimization_method", nargs="?", default="No") # args for bool
parser.add_argument("optimization_conf", nargs="?", default="No") # args for bool
parser.add_argument("optimization_lig", nargs="?", default="No") # args for bool
parser.add_argument("pre_optimization_lig", nargs="?", default="No") # args for bool
parser.add_argument("genconformer", nargs="?", default="No") # args for bool
parser.add_argument("nprocs", type=int, default=nprocs_all)
parser.add_argument("thr_fmax", type=float, default=0.05)
parser.add_argument("maxiter", type=float, default=500)

parser.add_argument("num_conformers", type=int, default=50)
parser.add_argument("max_attempts", type=int, default=100)
parser.add_argument("prune_rms_thresh", type=float, default=0.2)
parser.add_argument("opt_prune_rms_thresh", type=float, default=0.2)


def getBoolStr(string):
    string = string.lower()
    if "true" in string or "yes" in string:
        return True
    elif "false" in string or "no" in string:
        return False
    else:
        print("%s is bad input!!! Must be Yes/No or True/False" %string)
        sys.exit(1)


def setG16calculator(lig, file_base, label, WORK_DIR):
    lig.setG16Calculator(
            label="%s/g16_%s/%s"%(WORK_DIR, label, file_base),
            chk="%s.chk"%file_base,
            nprocs=nprocs,
            xc="HF",
            basis="6-31g*",
            scf="maxcycle=100",
    )
    return lig


def clusterConf(conf_dir, rmsd_thresh):
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
            Chem.rdMolAlign.AlignMol(mol1, mol2)
            dist_matrix[i, j] = Chem.rdMolAlign.CalcRMS(mol1, mol2)

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


def pruneConfs(cluster_conf, confs_energies, conf_dir):
    for fl_names in cluster_conf.values():
        for i, fl_name in enumerate(fl_names):
            e = float(confs_energies.loc[confs_energies["FileName"] == fl_name, " Energy(eV)"].item())
            if i == 0:
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


def runConfGen(file_name):
    "Starting ligand preparetion process... "
    mol_path= "%s/%s"%(structure_dir, file_name)

    file_base = file_name.split(".")[0]
    #create destination directory
    WORK_DIR = file_base
    if os.path.exists(WORK_DIR):
        shutil.rmtree(WORK_DIR)
    os.mkdir(WORK_DIR)

    #Flags
    # default mm calculator set to False
    mmCalculator=False
    # default adding H is False
    addH = False

    # if desire adding H by openbabel
    prefix = ""
    if add_hydrogen:
        addH = True
        prefix += "addH_"
    if optimization_lig or optimization_conf:
        prefix += "opt_"

    # initialize confGen
    lig = ligPrep(mol_path, addH, WORK_DIR)
    lig.setOptMethod(optimization_method)
    #  lig.writeRWMol2File("test/test.xyz")

    if "ani2x" in calculator_type.lower():
        lig.setANI2XCalculator()
    elif "g16" in calculator_type.lower():
        lig = setG16calculator(lig, file_base, label="calculation", WORK_DIR=WORK_DIR)
    elif "uff" in calculator_type.lower():
        if optimization_conf:
            print("UFF calculator not support optimization")
            sys.exit(1)
        else:
            mmCalculator=True

    # set optimizetion parameters
    lig.setOptParams(fmax=thr_fmax, maxiter=1000)

    if pre_optimization_lig:
        print("G16 Optimization process.. before generations")
        lig.geomOptimization()

    if genconformer:
        out_file_path="%s/%sminE_conformer.sdf"%(WORK_DIR, prefix)
        lig.genMinEGonformer(
            file_path=out_file_path,
            numConfs=num_conformers,
            maxAttempts=max_attempts,
            pruneRmsThresh=prune_rms_thresh,
            mmCalculator=mmCalculator,
            optimization_conf=optimization_conf,
        )

        if optimization_conf:
            CONF_DIR = WORK_DIR + "/conformers"
            confs_energies = pd.read_csv(f"{CONF_DIR}/confs_energies.csv")
            cluster_conf = clusterConf(CONF_DIR, rmsd_thresh=opt_prune_rms_thresh)
            if cluster_conf != 0:
                pruneConfs(cluster_conf, confs_energies, CONF_DIR)

        print("Conformer generation process is done")
        if not optimization_conf and optimization_lig:
            print("Optimization for minumum energy conformer")
            lig.geomOptimization()

    else:
        out_file_path="%s/%s%s.sdf"%(WORK_DIR, prefix, file_base)
        # geometry optimizaton for ligand
        if  optimization_lig:
            #  ase_atoms = lig.rwMol2AseAtoms()
            lig.geomOptimization()

    # write minimun energy conformer to sdf file
    lig.writeRWMol2File(out_file_path)

    #  for the bug of reading sfd file which have charges in ase
    try:
        atoms = read(out_file_path)
    except:
        out_file_path="%s/%s%s.xyz"%(WORK_DIR, prefix, file_base)
        lig.writeRWMol2File(out_file_path)
        atoms = read(out_file_path)


if __name__ == "__main__":
    args = parser.parse_args()
    structure_dir = args.structure_dir
    calculator_type = args.calculator_type

    optimization_method = args.optimization_method

    optimization_conf = getBoolStr(args.optimization_conf)
    optimization_lig = getBoolStr(args.optimization_lig)
    pre_optimization_lig = getBoolStr(args.pre_optimization_lig)
    genconformer = getBoolStr(args.genconformer)
    add_hydrogen = getBoolStr(args.add_hydrogen)

    nprocs = args.nprocs
    thr_fmax = args.thr_fmax
    maxiter = args.maxiter

    #get conformer generator parameters
    num_conformers = args.num_conformers
    max_attempts = args.max_attempts
    prune_rms_thresh = args.prune_rms_thresh
    opt_prune_rms_thresh = args.opt_prune_rms_thresh

    file_names = [item for item in os.listdir(structure_dir) if not item.startswith(".")]
    failed_csv = open("failed_files.csv", "w")
    failed_csv.write("FileNames,\n")

    for file_name in file_names:
        file_base = file_name.split(".")[0]
        #  try:
        runConfGen(file_name)
        #  except:
        #      print("Error for %s file !!! Skipping...")
        #      failed_csv.write(file_name+",\n")
        break
    failed_csv.close()

