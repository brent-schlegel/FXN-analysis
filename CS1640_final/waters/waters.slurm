#!/bin/bash
#SBATCH --job-name=ISCU_waters
#SBATCH -p dept_cpu
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36


python ../../md-scripts/waters.py ../../sim/complexwfe_M140I.prmtop ../dcd/complexwfe_M140I_1.dcd waters_1.csv
python ../../md-scripts/waters.py ../../sim/complexwfe_M140I.prmtop ../dcd/complexwfe_M140I_2.dcd waters_2.csv
python ../../md-scripts/waters.py ../../sim/complexwfe_M140I.prmtop ../dcd/complexwfe_M140I_3.dcd waters_3.csv
python ../../md-scripts/waters.py ../../sim/complexwfe_M140I.prmtop ../dcd/complexwfe_M140I_4.dcd waters_4.csv
python ../../md-scripts/waters.py ../../sim/complexwfe_M140I.prmtop ../dcd/complexwfe_M140I_5.dcd waters_5.csv
python ../../md-scripts/waters.py ../../sim/complexwfe_M140I.prmtop ../dcd/complexwfe_M140I_6.dcd waters_6.csv
python ../../md-scripts/waters.py ../../sim/complexwfe_M140I.prmtop ../dcd/complexwfe_M140I_7.dcd waters_7.csv
python ../../md-scripts/waters.py ../../sim/complexwfe_M140I.prmtop ../dcd/complexwfe_M140I_8.dcd waters_8.csv
python ../../md-scripts/waters.py ../../sim/complexwfe_M140I.prmtop ../dcd/complexwfe_M140I_9.dcd waters_9.csv
python ../../md-scripts/waters.py ../../sim/complexwfe_M140I.prmtop ../dcd/complexwfe_M140I_10.dcd waters_10.csv
