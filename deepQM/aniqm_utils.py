#
# Import ANI ensemble loader
import sys
#sys.path.append('/home/olexandr/notebooks/ASE_ANI/lib')
from ase_interface import ANIENS
from ase_interface import aniensloader

import  ase
import time
from ase import units
from ase.io import read, write, trajectory
from ase.optimize import BFGS, LBFGS

import pandas as pd
import os, re
from operator import itemgetter

import numpy  as np

def pdb2xyz(structure_path, file_base):
    cmd = "obabel %s/%s.pdb -O %s/%s.xyz 2> /dev/null" %(structure_path, file_base, structure_path, file_base)
    os.system("%s" %cmd)

def prepare_xyz_complex(work_path, file_base):
    fl_xyz = open("{}/{}.xyz".format(work_path, file_base), "w")
    #fl_xyz = open("{}.xyz".format(key_word), "w")
    data_xyz = []
    lines = open("{}/{}.pdb".format(work_path, file_base)).readlines()
    for line in lines:
        #print(line)
        if "ATOM" in line or "HETATM" in line:
            data_xyz.append(line)
    init_line = len(data_xyz)
    fl_xyz.write(str(init_line)+"\n{}/{}.xyz\n".format(work_path, file_base)) # file root to second row
    for line in data_xyz:
        split_line = line.split()
        line = "\t".join(itemgetter(5, 6, 7)(split_line)) # select data for xyz format and convert from tuple to str
        atom_sym = "".join(re.findall("[a-zA-Z]+", split_line[2]))[0] # extract atom symbol 
        line = atom_sym + "\t" + line

        fl_xyz.write(line)
        fl_xyz.write("\n")
    fl_xyz.close()

def prepare_xyz_grp(key_word, work_path, file_base):
    fl_xyz = open("{}/{}_{}.xyz".format(work_path, file_base, key_word), "w")
    #fl_xyz = open("{}.xyz".format(key_word), "w")
    data_xyz = []
    lines = open("{}/{}.pdb".format(work_path, file_base)).readlines()
    for line in lines[3:]:
        if key_word in line:
            data_xyz.append(line)
    init_line = len(data_xyz)
    fl_xyz.write(str(init_line)+"\n{}/{}_{}.xyz\n".format(work_path, file_base, key_word)) # file root to second row
    for line in data_xyz:
        split_line = line.split()
        line = "\t".join(itemgetter(5, 6, 7)(split_line)) # select data for xyz format and convert from tuple to str
        atom_sym = "".join(re.findall("[a-zA-Z]+", split_line[2]))[0] # extract atom symbol 
        line = atom_sym + "\t" + line

        fl_xyz.write(line)
        fl_xyz.write("\n")
    fl_xyz.close()

