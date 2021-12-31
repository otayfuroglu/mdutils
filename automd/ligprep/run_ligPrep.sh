#! /usr/bin/env bash


ligPrep_DIR="/cta/users/otayfuroglu/workspace/mdutils/automd/ligprep"
PYTHON_DIR="$HOME/miniconda3/bin"



struct_dir=test

# adding hydrogen if missing (yes/no)
add_hydrogen=yes

# generate conformer if desired (yes/no)
genconformer=yes

#configuration for conformer generator parameters
num_conformers=5
max_attempts=1000
prune_rms_thresh=0.2

# select caclulator type (ani2x/g16) for optimization conf
# caculator_type=g16
caculator_type="ani2x"

# perform geometry optimization for conformers if desired (yes/no)
optimization_conf=no

# perform geometry optimization for orginal ligand if desired (yes/no)
optimization_lig=yes

# set thrshold fmax for optimization (default=0.01)
thr_fmax=0.7

#maximum iteration for optimization
maxiter=500


$PYTHON_DIR/python $ligPrep_DIR/runLigPrep.py $struct_dir $add_hydrogen $caculator_type\
	$optimization_conf $optimization_lig $genconformer $thr_fmax $maxiter $num_conformers $max_attempts $prune_rms_thresh



