
from ligPrep import ligPrep
<<<<<<< HEAD
mol_path = "./test_lig.xyz"
=======
mol_path = "./test/M20.sdf"
>>>>>>> f999c4ac4c96399756b82a80c3b6524f6ba2df26
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
