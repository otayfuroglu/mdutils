
from ligPrep import ligPrep
mol_path = "./test_lig.xyz"
lig = ligPrep(mol_path)
#  lig.obMol2RWmol()
lig.writeRWMol2File("test_rw_file.pdb")
lig.genMinEGonformer("minE_conformer.sdf")




#  lig.convertFileFormat("mol2", "./test/convetedFile.mol2")
#  lig.writeOBMol2File("mol", "test/test_removeH.mol")
#  lig.addHWithOB()
#  lig.writeOBMol2File("mol", "test/test.mol")
#  lig.getObmolFromRWmol()
#  lig.writeRWMol2File("./test/test.pdb")
#  lig.printTest()
