#!/bin/sh

# AUTHOR: ALEX GIL GOMEZ
# INPUTS: TXT FILE WITH KEGG SPECIES NAMES. For example, bac_species.txt or fun_species.txt
# OUTPUT: SPECIES TREES
# INFO: GENERATES SPECIES TREES USING PROTEIN ALIGNMENTS AND 18/16S RNA GIVEN A SPECIES LIST. 
# HOW TO RUN: CHANGE VARIABLE Filename, AND RUN USING sbatch 0.PhySpeTree.sh
# DEPENDENCIES: REQUIRES THE PROGRAM PHYSPETREE INSTALL IN CONDA ENVIRONMENT. https://pypi.org/project/PhySpeTree/

###MODULES


#ENVIRONMENTS
conda activate base

# WORKING DIRECTORY

cd ../data/2.1.Phylogenetics_PhySpeTree

# VARIABLE. Filename (input txt file with KEGG species names). For example bac_species.txt or fun_species.txt

Filename="bac_sp.txt"

PhySpeTree autobuild -i bac_sp.txt -o bac_sp.txt.hcp --iqtree --hcp

PhySpeTree autobuild -i bac_sp.txt -o bac_sp.txt.srna --iqtree --srna

