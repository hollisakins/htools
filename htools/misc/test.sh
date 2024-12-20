#!/bin/bash
#SBATCH --job-name=cw_inspecs
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks=48
#SBATCH --ntasks-per-node=48
#SBATCH --time=20:00:00
#SBATCH --output=/home/hakins/cw_inspecs.%j

source /home/hakins/.bash_profile
conda activate base

python -m htools.misc.fitsmap_plots

exit 0