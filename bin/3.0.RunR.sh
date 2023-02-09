#!/bin/sh

#SBATCH --nodes=1
#SBATCH -p long-40core
#SBATCH --time=2-00:00:00

# AUTHOR: ALEX GIL GOMEZ
# INPUTS: Dataset DDI
# OUTPUT: list with tSNE results
# DEPENDENCIES: REQUIRES the R packages bigMap and bigmemory

###MODULES
module load R
#conda activate r_env

#ENVIRONMENTS

##RUN
#sbatch --job-name=jobname --ouutput=output --export=from=5,to=505 script.sh

# WORKING DIRECTORY. Should be ch3-netbio
pwd


#Rscript bin/3.1.tSNE_BigMap.R $from $to
Rscript bin/3.3.Modularity_test.R
#Rscript bin/3.4.Modularity_test_comparison.R
