
from ligPrep import ligPrep
import argparse
import os, sys



parser = argparse.ArgumentParser(description="Give something ...")
parser.add_argument("structure_dir", type=str)
parser.add_argument("add_hydrogen", nargs="?", default=False) # args for bool
parser.add_argument("calculator_type", type=str)
parser.add_argument("optimization", nargs="?", default=False) # args for bool
parser.add_argument("thr_fmax", type=float, default=0.05)
parser.add_argument("maxiter", type=float, default=500)

parser.add_argument("num_conformers", type=int, default=50)
parser.add_argument("max_attempts", type=int, default=100)
parser.add_argument("prune_rms_thresh", type=float, default=0.2)

args = parser.parse_args()
structure_dir = args.structure_dir
calculator_type = args.calculator_type

optimization = args.optimization.lower()
if "true" in optimization or "yes" in optimization:
    optimization = True
else:
    optimization = False

add_hydrogen = args.add_hydrogen.lower()
if "true" in add_hydrogen or "yes" in add_hydrogen:
    add_hydrogen = True
else:
    add_hydrogen = False

thr_fmax = args.thr_fmax
maxiter = args.maxiter


#get conformer generator parameters
num_conformers = args.num_conformers
max_attempts = args.max_attempts
prune_rms_thresh = args.prune_rms_thresh

def setG16calculator(lig, file_base):
    lig.setG16Calculator(
            label="tmp/g16_%s"%file_base,
            chk="g16_%s.chk"%file_base,
            xc="B3LYP",
            basis="sto-3g",
            scf="maxcycle=100",
            multiplicity=1,
            extra="Pop=(MK)",
    )
    return lig


def run(mol_path):
    file_base = mol_path.split("/")[-1].replace(".mol2", "")

    # initialize ligPrep
    lig = ligPrep(mol_path)
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
        if optimization:
            print("UFF calculator not support optimization")
            sys.exit(1)
        else:
            mmCalculator=True

    if optimization:
        lig.setOptParams(fmax=0.05, maxiter=1000)

    out_file_path="minE_conformer_%s.sdf"%file_base
    lig.genMinEGonformer(file_path=out_file_path,
                         numConfs=num_conformers,
                         maxAttempts=max_attempts,
                         pruneRmsThresh=prune_rms_thresh,
                         mmCalculator=mmCalculator,
                         optimization=optimization,
                        )

    if sp4esp:
        from ase.io import read
        atoms = read(out_file_path)

        lig = setG16calculator(lig, file_base)
        lig.calcSPEnergy(atoms)

file_names = os.listdir(structure_dir)
for file_name in file_names:
    mol_path= "%s/%s"%(structure_dir, file_name)
    run(mol_path)
#
#      start = time.time()
#
#      mol_path= "./test/untested/%s"%file_name
#      lig = ligPrep(mol_path)
#      lig.setMaxCycle(250)
#      lig.genMinEGonformer("minE_conformer/minE_conformer_%s"%file_name)
#
#      spend_time =(time.time() - start) / 60.0
#      print("time: {0:.2f} minute".format(spend_time))
#

#  lig.convertFileFormat("mol2", "./test/convetedFile.mol2")
#  lig.writeOBMol2File("mol", "test/test_removeH.mol")
#  lig.addHWithOB()
#  lig.writeOBMol2File("mol", "test/test.mol")
#  lig.getObmolFromRWmol()
#  lig.writeRWMol2File("./test/test.pdb")
#  lig.printTest()
