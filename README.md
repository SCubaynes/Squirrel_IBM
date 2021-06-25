# Squirrel_IBM
Script and simulations to run IBMs and calculate kinship structure for the ground squirrel example.
Results are published in "The Kinship Matrix: Inferring the Kinship Structure of a Population From Its Demography", published in 2021 in Ecology Letters from authors Christophe F. D. Coste, François Bienvenu, Victor Ronget, Sarah Cubaynes, Juan Pablo Ramirez Loza and Samuel Pavard.

Results are saved in R.Data files :
- kinship, a list containing kinship matrices for each replicate
- lambda, a containing growth rate at last iteration for each replicate
- w, a list with stable age structure calculated at the last iteration for each replicate 
- popSimRep,an list of agentmatrix objects containing the simulated populations - require package NetLogoR
- popSimRep_matrixformat, a list of matrix objects containing the simulated populations - does not require package NetLogoR
