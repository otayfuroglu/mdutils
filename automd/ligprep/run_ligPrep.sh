#! /usr/bin/env bash


ligPrep_DIR="/cta/users/otayfuroglu/workspace/mdutils/automd"
PYTHON_DIR="$HOME/miniconda3/bin"



struct_dir=test0

# adding hydrogen if missing (yes/no)
add_hydrogen=yes

# select caclulator type (ani2x/g16)
# caculator_type=g16
caculator_type=ani2x

# performe geometry optimization if desired (yes/no)
optimization=no


# set thrshold fmax for optimization (default=0.01)
thr_fmax=1000.01

#maximum iteration for optimization
maxiter=500

#configuration for conformer generator parameters
num_conformers=50
max_attempts=1000
prune_rms_thresh=0.2

$PYTHON_DIR/python $ligPrep_DIR/runLigPrep.py $struct_dir $add_hydrogen $caculator_type $optimization $thr_fmax $maxiter $num_conformers $max_attempts $prune_rms_thresh



