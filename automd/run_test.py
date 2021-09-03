
from ligPrep import ligPrep
import time
import os

file_names = os.listdir("./test/untested")


for file_name in file_names:

    start = time.time()

    mol_path= "./test/untested/%s"%file_name
    lig = ligPrep(mol_path)
    lig.setMaxCycle(250)
    lig.genMinEGonformer("minE_conformer/minE_conformer_%s"%file_name)

    spend_time =(time.time() - start) / 60.0
    print("time: {0:.2f} minute".format(spend_time))



#  lig.convertFileFormat("mol2", "./test/convetedFile.mol2")
#  lig.writeOBMol2File("mol", "test/test_removeH.mol")
#  lig.addHWithOB()
#  lig.writeOBMol2File("mol", "test/test.mol")
#  lig.getObmolFromRWmol()
#  lig.writeRWMol2File("./test/test.pdb")
#  lig.printTest()
