
from ligPrep import ligPrep
import argparse
import os, sys



parser = argparse.ArgumentParser(description="Give something ...")
parser.add_argument("structure_dir", type=str)
parser.add_argument("add_hydrogen", nargs="?", default="No") # args for bool
parser.add_argument("calculator_type", type=str)
parser.add_argument("optimization_conf", nargs="?", default="No") # args for bool
parser.add_argument("optimization_lig", nargs="?", default="No") # args for bool
parser.add_argument("genconformer", nargs="?", default="No") # args for bool
parser.add_argument("thr_fmax", type=float, default=0.05)
parser.add_argument("maxiter", type=float, default=500)

parser.add_argument("num_conformers", type=int, default=50)
parser.add_argument("max_attempts", type=int, default=100)
parser.add_argument("prune_rms_thresh", type=float, default=0.2)

args = parser.parse_args()
structure_dir = args.structure_dir
calculator_type = args.calculator_type

def getBoolStr(string):
    string = string.lower()
    if "true" in string or "yes" in string:
        return True
    elif "false" in string or "no" in string:
        return False
    else:
        print("%s is bad input!!! Must be Yes/No or True/False" %string)
        sys.exit(1)


optimization_conf = getBoolStr(args.optimization_conf)
optimization_lig = getBoolStr(args.optimization_lig)
genconformer = getBoolStr(args.genconformer)
add_hydrogen = getBoolStr(args.add_hydrogen)

thr_fmax = args.thr_fmax
maxiter = args.maxiter


#get conformer generator parameters
num_conformers = args.num_conformers
max_attempts = args.max_attempts
prune_rms_thresh = args.prune_rms_thresh

def setG16calculator(lig, file_base, WORK_DIR):
    lig.setG16Calculator(
            label="%s/g16_calculation/%s"%(WORK_DIR, file_base),
            chk="%s.chk"%file_base,
            xc="HF",
            basis="sto-3g",
            scf="maxcycle=100",
            extra="Pop=(MK) IOP(6/50=1)",
            addsec="%s.esp"%file_base,
    )
    return lig


def runLigPrep(file_base):
    "Starting ligand preparetion process... "
    mol_path= "%s/%s.mol2"%(structure_dir, file_base)

    #create destination directory
    WORK_DIR = file_base
    if not os.path.exists(WORK_DIR):
        os.mkdir(WORK_DIR)

    # initialize ligPrep
    lig = ligPrep(mol_path, WORK_DIR)
    #  lig.writeRWMol2File("test/test.xyz")

    # if desire adding H by rdkit
    if add_hydrogen:
        lig.addHwithRD()

    # defaul mm calculator set to False
    mmCalculator=False

    sp4esp = True
    if "ani2x" in calculator_type.lower():
        lig.setANI2XCalculator()
    elif "g16" in calculator_type.lower():
        sp4esp = False
        lig = setG16calculator(lig, file_base)
    elif "uff" in calculator_type.lower():
        if optimization_conf:
            print("UFF calculator not support optimization")
            sys.exit(1)
        else:
            mmCalculator=True

    #  if optimization_conf:
    # set optimizetion parameters
    lig.setOptParams(fmax=thr_fmax, maxiter=1000)


    prefix = ""
    if optimization_lig or optimization_conf:
        prefix = "opt_"
    out_file_path="%s/%sminE_conformer.xyz"%(WORK_DIR, prefix)
    if genconformer:
        ase_atoms = lig.genMinEGonformer(
            file_path=out_file_path,
            numConfs=num_conformers,
            maxAttempts=max_attempts,
            pruneRmsThresh=prune_rms_thresh,
            mmCalculator=mmCalculator,
            optimization_conf=optimization_conf,
        )

        print("Conformer generation process is done")
        print("Selected minimun energy conformer")
        if  not optimization_conf and optimization_lig:
            print("Optimization for minumum energy conformer")
            lig.geomOptimization(ase_atoms)
            lig.writeAseAtoms(out_file_path)
        else:
            lig.writeAseAtoms(out_file_path)

    else:
        # geometry optimizaton for ligand
        if  optimization_lig:
            ase_atoms = lig.rwMol2AseAtoms()
            lig.geomOptimization(ase_atoms)
        # ligad coordintes write to xyz file
        lig.writeAseAtoms(out_file_path)


    if sp4esp:
        from ase.io import read
        atoms = read(out_file_path)

        lig = setG16calculator(lig, file_base, WORK_DIR)
        lig.calcSPEnergy(atoms)

def main():
    file_names = [item for item in os.listdir(structure_dir) if item.endswith(".mol2")]
    for file_name in file_names:
        file_base = file_name.replace(".mol2", "")
        runLigPrep(file_base)
        os.system("bash ligPrep_utils.sh %s" %file_base)

if __name__ == "__main__":
    main()
