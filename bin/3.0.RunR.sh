#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=extended-96core
#SBATCH --time=7-00:00:00

# AUTHOR: ALEX GIL GOMEZ
# INPUTS: Dataset DDI
# OUTPUT: list with tSNE results
# DEPENDENCIES: REQUIRES the R packages bigMap and bigmemory

###MODULES
module load R
#conda activate r_env

#ENVIRONMENTS

##RUN
#sbatch --job-name=jobname --output=output --export=from=5,to=505 script.sh

# WORKING DIRECTORY. Should be ch3-netbio
pwd


#Rscript bin/3.1.tSNE_BigMap.R $from $to
# Rscript bin/3.3.Modularity_test.R allddi ppx_205 ppx_1455
# Rscript bin/3.3.Modularity_test.R sen2in1 ppx_115 ppx_455
Rscript bin/3.4.Modularity_test_comparison.R
