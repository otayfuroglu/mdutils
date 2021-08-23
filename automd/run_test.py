
from ligPrep import ligPrep
mol_path = "./test/M20.sdf"
lig = ligPrep(mol_path)
lig.genMinEGonformer("minE_conformer.sdf")




#  lig.convertFileFormat("mol2", "./test/convetedFile.mol2")
#  lig.writeOBMol2File("mol", "test/test_removeH.mol")
#  lig.addHWithOB()
#  lig.writeOBMol2File("mol", "test/test.mol")
#  lig.getObmolFromRWmol()
#  lig.writeRWMol2File("./test/test.pdb")
#  lig.printTest()
