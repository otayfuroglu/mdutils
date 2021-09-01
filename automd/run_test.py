
from ligPrep import ligPrep
import time

start = time.time()

mol_path= "./test/M20.sdf"
lig = ligPrep(mol_path)
lig.genMinEGonformer("minE_conformer_M20.sdf")

spend_time =(time.time() - start) / 60.0
print("time: {0:.2f} minute".format(spend_time))



#  lig.convertFileFormat("mol2", "./test/convetedFile.mol2")
#  lig.writeOBMol2File("mol", "test/test_removeH.mol")
#  lig.addHWithOB()
#  lig.writeOBMol2File("mol", "test/test.mol")
#  lig.getObmolFromRWmol()
#  lig.writeRWMol2File("./test/test.pdb")
#  lig.printTest()
