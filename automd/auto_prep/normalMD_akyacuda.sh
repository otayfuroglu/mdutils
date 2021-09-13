source /truba/sw/centos7.3/comp/intel/PS2018-update2/bin/compilervars.sh intel64
module load centos7.3/app/gromacs/2020-impi-mkl-PS2018-GOLD-CUDA
module load centos7.3/comp/cmake/3.10.1
module load centos7.3/comp/gcc/6.4
module load centos7.3/lib/cuda/10.0

#export GMX_GPU_DD_COMMS=true
#export GMX_GPU_PME_PP_COMMS=true
#export GMX_FORCE_UPDATE_DEFAULT_GPU=true

FREE_ENERGY=$(pwd)
MDP=$FREE_ENERGY/MDP
MDRUNS=$FREE_ENERGY/MDRUNS
export GMX_MAXBACKUP=-1
export GMX_MAXWARN=-1
output=$$.out
exec>$output 2>&1

LAMBDA=0

cd $FREE_ENERGY
#################################
# ENERGY MINIMIZATION 1: STEEP  #
#################################
echo "Starting minimization for lambda = $LAMBDA..." 
cd $MDRUNS/EM
# Iterative calls to grompp and mdrun to run the simulations
gmx grompp -f $MDP/EM/em_steep.mdp -c $FREE_ENERGY/solv_ions.gro -r $FREE_ENERGY/solv_ions.gro -p $FREE_ENERGY/topol.top -o min_steep_$LAMBDA.tpr -n $FREE_ENERGY/index.ndx 
gmx mdrun -v -deffnm min_steep_$LAMBDA 
sleep 3

for LAMBDA in {0..2..1};do


#######################
## NVT EQUILIBRATION #
######################
echo "Starting constant volume equilibration..."
cd $MDRUNS/NVT
gmx grompp -f $MDP/NVT/nvt.mdp -c $MDRUNS/EM/min_steep_0.gro -r $MDRUNS/EM/min_steep_0.gro -p $FREE_ENERGY/topol.top -o nvt_$LAMBDA.tpr -n $FREE_ENERGY/index.ndx
gmx mdrun -v -ntmpi 4 -ntomp 10 -npme 1 -nb gpu -bonded gpu -pme gpu -deffnm nvt_$LAMBDA 
echo "Constant volume equilibration complete."
sleep 3
#
######################
## NPT part1 EQUILIBRATION #
######################
echo "Starting constant pressure equilibration using Berendsen barostat..."
cd $MDRUNS/NPT
gmx grompp -f $MDP/NPT/npt.mdp -c $MDRUNS/NVT/nvt_$LAMBDA.gro -r $MDRUNS/NVT/nvt_$LAMBDA.gro -p $FREE_ENERGY/topol.top -t $MDRUNS/NVT/nvt_$LAMBDA.trr -o npt_$LAMBDA.tpr -n $FREE_ENERGY/index.ndx 
gmx mdrun  -v -ntmpi 4 -ntomp 10 -npme 1 -nb gpu -bonded gpu -pme gpu -deffnm npt_$LAMBDA 
echo "Constant pressure equilibration complete with Berendsen barostat."
sleep 10

######################
### NPT part2 EQUILIBRATION #
#######################
echo "Starting constant pressure equilibration using Parrinello-Rahman barostat..."
cd $MDRUNS/NPT
gmx grompp -f $MDP/NPT/npt2.mdp -c $MDRUNS/NPT/npt_$LAMBDA.gro -r $MDRUNS/NPT/npt_$LAMBDA.gro -p $FREE_ENERGY/topol.top -t $MDRUNS/NPT/npt_$LAMBDA.trr -o npt2_$LAMBDA.tpr -n $FREE_ENERGY/index.ndx 
gmx mdrun  -v -ntmpi 4 -ntomp 10 -npme 1 -nb gpu -bonded gpu -pme gpu -deffnm npt2_$LAMBDA 
echo "Constant pressure equilibration complete with Parrinello-Rahman barostat."
sleep 3 

####################
#### PRODUCTION MD #
####################
echo "Starting production MD simulation..."
cd $MDRUNS/Production_MD
gmx grompp -f $MDP/Production_MD/md.mdp -c $MDRUNS/NPT/npt2_$LAMBDA.gro -r $MDRUNS/NPT/npt2_$LAMBDA.gro  -p $FREE_ENERGY/topol.top -t $MDRUNS/NPT/npt2_$LAMBDA.cpt -o md_$LAMBDA.tpr -n $FREE_ENERGY/index.ndx 
gmx mdrun -v -ntmpi 4 -ntomp 10 -npme 1 -nb gpu -bonded gpu -pme gpu -deffnm md_$LAMBDA 
echo "Production MD complete."
echo "Ending. Job completed for lambda = $LAMBDA"
sleep 3
done

