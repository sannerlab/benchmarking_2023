![stripSet](https://github.com/sannerlab/benchmarking_2023/blob/main/Fnat2.png)

This folder contains python scripts and required data to perform calculations for the paper "Predicting protein-peptide interactions: benchmarking Deep Learning techniques and a comparison with focused docking ".
Please cite: 


The function "estimate_energies_for_pdb" from OpenMM_functions.py can be used perform minimization and various energy calculations for a protein-peptide complex.

The percentage docking success rate as a function of Fnat cutoff can be calculated by "plot_fnat_vs_success_shared.py".  

To compare the performance of different methods "cross_performance_analysis.py" script can be used.

The details of PDB files, protein sequences, and peptide sequences for benchmarking is given in  "pdb_list.txt" and "pdb_sequences.csv" files.

The file "all_fnat.dat" contains Fnat values from all methods for different top solutions. 
