#!/bin/sh

#SBATCH --job-name=tSNE
#SBATCH -o tSNE.sh.out
#SBATCH --nodes=1
#SBATCH -p extended-40core
#SBATCH --time=7-00:00:00

# AUTHOR: ALEX GIL GOMEZ
# INPUTS: Dataset DDI
# OUTPUT: list with tSNE results
# INFO: Runs the script for 4.2.tSNE  
# HOW TO RUN: sbatch nameofscript.sh
# DEPENDENCIES: REQUIRES the R packages bigMap and bigmemory

###MODULES


#ENVIRONMENTS


# WORKING DIRECTORY. Should be ch3-netbio
pwd


# VARIABLE. 
Rscript code/4.2.tSNE_BigMap.R

