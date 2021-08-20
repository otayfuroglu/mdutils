
from ligPrep import ligPrep
mol_path = "./Test_sdf_files/M20.sdfz"
lig = ligPrep(mol_path)
lig.addH()
lig.obMol2RWmol()
lig.writeRWMol2File("test_rw_file.mol2")
#  lig.removeH()
#  lig.addH()
#  lig.convertFileFormat("mol2", "./test/convetedFile.mol2")
#  lig.writeOBMol2File("mol", "test/test_removeH.mol")
#  lig.addHWithOB()
#  lig.writeOBMol2File("mol", "test/test.mol")
#  lig.getObmolFromRWmol()
#  lig.writeRWMol2File("./test/test.pdb")
#  lig.printTest()
