#! /usr/bin/python3

from ase import Atoms
from ase.io import write
import os, io
import numpy as np
import shutil
import argparse



def main(xyzDIR, outFileName):
    #output file name

    fileNames = [item for item in os.listdir(xyzDIR) if item.endswith(".xyz")]
    print("Nuber of xyx files --> ", len(fileNames))

    for i, fileName in enumerate(fileNames):

        if i == 0:
            assert not os.path.exists(outFileName), "%s exists" %(outFileName)
        inFile = open(os.path.join(xyzDIR, fileName), "r")
        fileStr = inFile.read()
        outFile = open(outFileName, 'a')
        outFile.write(fileStr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-xyzDIR", "--xyzDIR", type=str, required=True, help="give xyz files directory")
    parser.add_argument("-o", "--outFileName", type=str, required=True, help="give out file name with extention")

    args = parser.parse_args()
    xyzDIR = args.xyzDIR
    outFileName = args.outFileName
    main(xyzDIR, outFileName)
