#!/bin/bash
#SBATCH --job-name=ISCU_contacts
#SBATCH -p dept_cpu
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32

python contacts.py MDAnalysis_supported_topolgy MDAnalysis_supported_trajectory csv_output_name.csv

python ../../md-scripts/contacts.py ../../sim/complexwfe_M140I.prmtop ../dcd/complexwfe_M140I_1.dcd contacts_1.csv
python ../../md-scripts/contacts.py ../../sim/complexwfe_M140I.prmtop ../dcd/complexwfe_M140I_2.dcd contacts_2.csv
python ../../md-scripts/contacts.py ../../sim/complexwfe_M140I.prmtop ../dcd/complexwfe_M140I_3.dcd contacts_3.csv
python ../../md-scripts/contacts.py ../../sim/complexwfe_M140I.prmtop ../dcd/complexwfe_M140I_4.dcd contacts_4.csv
python ../../md-scripts/contacts.py ../../sim/complexwfe_M140I.prmtop ../dcd/complexwfe_M140I_5.dcd contacts_5.csv
python ../../md-scripts/contacts.py ../../sim/complexwfe_M140I.prmtop ../dcd/complexwfe_M140I_6.dcd contacts_6.csv
python ../../md-scripts/contacts.py ../../sim/complexwfe_M140I.prmtop ../dcd/complexwfe_M140I_7.dcd contacts_7.csv
python ../../md-scripts/contacts.py ../../sim/complexwfe_M140I.prmtop ../dcd/complexwfe_M140I_8.dcd contacts_8.csv
python ../../md-scripts/contacts.py ../../sim/complexwfe_M140I.prmtop ../dcd/complexwfe_M140I_9.dcd contacts_9.csv
python ../../md-scripts/contacts.py ../../sim/complexwfe_M140I.prmtop ../dcd/complexwfe_M140I_10.dcd contacts_10.csv
