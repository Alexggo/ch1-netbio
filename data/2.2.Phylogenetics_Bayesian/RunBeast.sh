#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=long-24core
#SBATCH --time=2-00:00:00

#RUN WITH:
#sbatch --export=XML_file='DDI_BD_rce' --job-name=${XML_file}.r --output=${XML_file}.out.txt RunBeast.sh

### LOAD MODULES
module load beast/2.6.7
cd ${XML_file}
beast ${XML_file}.xml
