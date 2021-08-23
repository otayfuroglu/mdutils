
from ligPrep import ligPrep
import time

start = time.time()

mol_pathc= "./test/M20.sdf"
lig = licPrep(mol_path)
lig.genMcnEGonformer("minE_conformer.sdf")

spend_time =(time.time() - start) / 60.0
print("time: {0:.2f} minute".format(spend_time))



#  lig.convertFileFormat("mol2", "./test/convetedFile.mol2")
#  lig.writeOBMol2File("mol", "test/test_removeH.mol")
#  lig.addHWithOB()
#  lig.writeOBMol2File("mol", "test/test.mol")
#  lig.getObmolFromRWmol()
#  lig.writeRWMol2File("./test/test.pdb")
#  lig.printTest()
