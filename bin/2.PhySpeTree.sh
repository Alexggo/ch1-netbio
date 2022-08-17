#!/bin/sh

# AUTHOR: ALEX GIL GOMEZ
# INPUTS: TXT FILE WITH KEGG SPECIES NAMES. For example, bac_species.txt or fun_species.txt
# OUTPUT: SPECIES TREES
# INFO: GENERATES SPECIES TREES USING PROTEIN ALIGNMENTS AND 18/16S RNA GIVEN A SPECIES LIST. 
# HOW TO RUN: CHANGE VARIABLE Filename, AND RUN USING sbatch 0.PhySpeTree.sh
# DEPENDENCIES: REQUIRES THE PROGRAM PHYSPETREE INSTALL IN CONDA ENVIRONMENT. https://pypi.org/project/PhySpeTree/

###MODULES


#ENVIRONMENTS
conda activate bioinformatics_py3

# WORKING DIRECTORY

cd ../data/2.1.Phylogenetics_PhySpeTree

# VARIABLE. Filename (input txt file with KEGG species names). For example bac_species.txt or fun_species.txt

Filename=(bac_species.txt fun_species.txt)

for m in "${Filename[@]}"
do
echo "Analysis of $Filename"
PhySpeTree autobuild -i $m -o $m.hcp --iqtree  --hcp 

PhySpeTree autobuild -i $m -o $m.srna --iqtree --srna 
done
echo "All done!"
