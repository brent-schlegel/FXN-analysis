#!/bin/bash
#SBATCH --job-name=ISCU_rmsfs
#SBATCH -p dept_cpu
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32

python ../../md-scripts/rmsf.py ../../sim/complexwfe_M140I.prmtop ../dcd/complexwfe_M140I_ 10
