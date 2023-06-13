#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=extended-40core
#SBATCH --time=7-00:00:00

#RUN WITH:
#sbatch --export=name='Phylo' 2.1.RunBeast.sl

#RUN WITH:

### LOAD MODULES
cd /gpfs/projects/RestGroup/agilgomez/projects/ch3-netbio/data/2.2.Phylogenetics_Bayesian

# Resume interrupted job:
/gpfs/projects/RestGroup/agilgomez/tools/beast.2.7.4/bin/beast -resume -statefile ${name}.xml.state  ${name}.xml

#/gpfs/projects/RestGroup/agilgomez/tools/beast.2.7.4/bin/beast ${name}.xml