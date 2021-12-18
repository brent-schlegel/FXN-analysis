import MDAnalysis as md
import numpy as np
import pandas as pd
import py3Dmol
import sys

pdb = md.Universe("complexwfe_M140I.pdb")
pdb.add_TopologyAttr('tempfactors')
protein = pdb.select_atoms('protein')

csv = pd.read_csv(sys.argv[1])

for index in range(len(csv['resid'])):
	res = protein.residues[index]
	if csv.loc[index,'p_score'] < 0.005:
		protein.residues[index].atoms.tempfactors = csv.loc[index,'wt_mean']-csv.loc[index,'mutant_mean']

protein.write(sys.argv[2])
