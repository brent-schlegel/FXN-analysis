#Frataxin Analysis

!All requisite code is found in the CS1640_final directory, GH wouldn't let me delete the other ones for some reason.

This project was completed pursuant to the University of Pittsburgh CS1640 : Bioinformatics Software Design course, in conjunction with Dr. David R. Koes and members of his lab at Pitt. The primary goal was to quantify the molecular interactions of mutant complexes assosciated with the development and potential treatment of Friedriech's Ataxia, a metabolic disorder hindering the formation of iron-sulfur clusters in mitochondria. This was achieved through Molecular Dynamics simulation and visualization using PyMol, SLURM, pMEMD, Amber, CHARMm, and Anaconda, as well as statistical extrapolation and analysis with Python3, cpptraj, Jupyter Notebook, and MDAnalysis. Analysis focused on three distinct categories of molecular dynamics: Protein-protein contacts, protein-water contacts, and root mean-square fluctuation (RMSF). Ten simulations per 5 mutants (50 total) were ran on the Pitt School of Medicine's computational cluster, as an array of simulations organized with SLURM. Analysis was performed through remote access. 

##sim:
Main SLURM file used to run the simulations on the cluster. Jobs were ran simultaneously in an array. 

##waters, contacts, rmsfs:

      *Statisical extrapolation of simulations, concatenation and scoring done with Python3 (pandas,numpy), 2-sided student's t-test assuming equal variance for each        residue, performed between the wild type (unmutated) and mutant simulated data across 10 simulations. 
      *Significance values are added in a new column of the dataframe. 
      *Output in CSV format. 

##final_distributions:
        
      Sample output of Jupyter visualizations including probability density, mean trend variance, difference btwn wild type & mutant, mapping of significant residues

##compare:
      
      Generating the PDB format of the mutated complexes, by comparing the mutant and wild-type values for a specific residue, setting the B-factor for significant        (p<0.005>) residues to the value of the difference of the wild-type and mutant. 
      If not significant, B-factor values are 0. This allows for visualization using Py3Dmol and PyMol. 

##jupyter_notebooks (ipynb):

      *Statistical analysis and visualization of simulated data vs wild type, graphing with matplotlib. 
      *Graphs are color-coded for either significance or wt-mutant differentiation, according to the legends. 
      *Plots proceed in complexity for a given subset of data (waters, for example). Starting with probability density distributions and ultimately ending up in the          analysis of significance. 
      *Complex 3D modelling performed with Py3Dmol

    *Distributions: Primary statisitical analysis as described above

    *Complex Analysis: (MDAnalysis, pc) Plotting of RMSF values across the whole complex, looking for variation across the peaks to indicate a significant variation in values. PCA generated by mapping the principal components from the wild-type to the mutated sequence, as they have the same sample size and variance. This analysis was more experimental, the results were inconclusive and served as an exercise in more laborious computational methods.

    *trajectory_analysis: Code provided by our Group Lead to serve as a framework for the general workflow of the statistical analysis.



##NOT INCLUDED:
*Simulation output (very large files)
*Proprietary code provided by group lead, used for AMBER preprocessing of wild-type data for Molecular Dynamics simulation & standardization of simulation output

##Authors:

*Brent Schlegel (Dietrich School of Arts and Sciences, University of Pittsburgh)
    github: brent-schlegel

*Aaron Fairchild (Dietrich School of Arts and Sciences, University of Pittsburgh)

*Tao Sheng (School of Computing and Information, University of Pittsburgh)

*Isaac Kim (Dietrich School of Arts and Sciences, University of Pittsburgh)

*Angela Tseng (School of Computing and Information, University of Pittsburgh)


##Acknowledgments:

*Dr. David Koes (Department of Computational and Systems Biology, University of Pittsburgh)

*Lucas Morley (Dietrich School of Arts and Sciences, University of Pittsburgh)
