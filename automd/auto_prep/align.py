
import rdkit 
import rdkit.Chem
import rdkit.Chem.rdMolAlign
import argparse

def alignTwoMol(mol1, mol2):    
    return rdkit.Chem.rdMolAlign.AlignMol(mol1, mol2)


def mainAlingTwoMol():
      parser = argparse.ArgumentParser(description="Give something ...")
      parser.add_argument("-i", "--mol_path", type=str, required=True, help="")
      parser.add_argument("-r", "--mol_ref_path", type=str, required=True, help="")
      parser.add_argument("-o", "--output", type=str, required=True, help="")

      args=parser.parse_args()
      mol_path=args.mol_path
      mol_ref_path=args.mol_ref_path
      output=args.output

      mol = rdkit.Chem.rdmolfiles.MolFromPDBFile(mol_path,removeHs=False)
      mol_ref = rdkit.Chem.rdmolfiles.MolFromPDBFile(mol_ref_path,removeHs=False)
      rmsd = alignTwoMol(mol, mol_ref)
      print("RMSD-->", rmsd)
 
      rdkit.Chem.rdmolfiles.MolToPDBFile(mol, output)





mainAlingTwoMol()
