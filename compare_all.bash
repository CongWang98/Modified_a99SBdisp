#!/bin/bash   
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --output=R-%x.%j.out
#SBATCH --time=96:00:00  

gmxPath="/home/gridsan/congwang/software/gromacs/bin/gmx"
gromacsTopDir=$PWD

python compare_energy.py folder=01_water gmxPath=$gmxPath gromacsTopDir=$gromacsTopDir > 01_water_output.txt
python compare_energy.py folder=02_dhfr gmxPath=$gmxPath gromacsTopDir=$gromacsTopDir > 02_dhfr_output.txt
python compare_energy.py folder=03_dhfr_in_water gmxPath=$gmxPath gromacsTopDir=$gromacsTopDir > 03_dhfr_in_water_output.txt

rm \#*
rm *.edr
rm *.gro
rm *.log
rm *.tpr
rm *.trr
rm mdout.mdp

