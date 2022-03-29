#! /usr/bin/env bash

ligPrep_DIR="$HOME/Desktop/mdutils/automd/ligprep"
PYTHON_DIR="$HOME/miniconda3/bin"

struct_dir=test

# adding hydrogen if missing (yes/no)
add_hydrogen=no

# set optizetions methods whichs availble in ASE (BFGS, LBFGS, GPMin, FIRE, Berny)
optimization_method=fire

# optimization ligand if desired before conformer generation (yes/no)
pre_optimization_lig=no

# generate conformer if desired (yes/no)
genconformer=yes

#configuration for conformer generator parameters
num_conformers=4
max_attempts=5000
prune_rms_thresh=0.05
opt_prune_rms_thresh=1.5

# select caclulator type (ani2x/g16) for optimization conf
# caculator_type=g16
caculator_type="ani2x"

# perform geometry optimization for conformers if desired (yes/no)
optimization_conf=yes

# perform geometry optimization for orginal ligand if desired (yes/no)
optimization_lig=no

# set number of procssors for g16 calcultor (default=all cpu)
nprocs=32

# set thrshold fmax for optimization (default=0.01)
thr_fmax=0.2

#maximum iteration for optimization
maxiter=500


$PYTHON_DIR/python $ligPrep_DIR/runConfGen.py $struct_dir $add_hydrogen $caculator_type\
	$optimization_method $optimization_conf $optimization_lig $pre_optimization_lig $genconformer\
       	$nprocs $thr_fmax $maxiter $num_conformers $max_attempts $prune_rms_thresh $opt_prune_rms_thresh


