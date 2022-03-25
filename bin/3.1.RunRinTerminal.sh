#!/bin/sh

#SBATCH --job-name=tSNE
#SBATCH -o tSNE2.sh.out
#SBATCH --nodes=1
#SBATCH -p long-40core
#SBATCH --time=2-00:00:00

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
# Rscript bin/3.1.tSNE_BigMap.R
Rscript bin/5.1.Networks.R