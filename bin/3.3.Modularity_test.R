#!/usr/bin/env Rscript

#' ---
#' title: "3.2.Modularity_test"
#' author: "Alex Gil"
#' date: "1/26/2021"
#' output: html_document
#' 
#' 
#' 
## -------------------------------------------------------------------------
library(pacman)
p_load(geomorph,tidyverse,broom,OUwie,phytools,ape,ggx,geomorph)


full_dataset <- read_csv("data/3.InteractionScores_tSNE/all_ppx_large.csv")

## Testing modularity of clustering solutions.

# Read phylogenetic tree
species <- c("Escherichia_coli_K-12_ebw",
             "Escherichia_coli_O8_IAI1_ecr",
             "Salmonella_enterica_serovar_Typhimurium_LT2_stm",
             "Salmonella_enterica_serovar_Typhimurium_14028S_seo",
             "Pseudomonas_aeruginosa_PAO1_pae",
             "Pseudomonas_aeruginosa_UCBPP-PA14_pau")

treetype <- "bac_species.txt.hcp"
treefile <- file.path('data/2.2.Phylogenetics_Bayesian/DDI_BD_str',"concatenate.trees.tre")
tree_nw <- read.nexus(treefile)
# Drop the species that are not needed.
not.species <- tree_nw$tip.label[!(tree_nw$tip.label %in% species)]
tree_nw <- drop.tip(tree_nw,not.species)

tree_nw$tip.label <- c("ebw","ecr","pae","pau","seo","stm")

mat <- full_dataset[,19:24]%>% as.matrix() |> t()


#test ppx 205-2255
mod_test <- list()
1:205
for (i in 1:2){
  labels <- full_dataset[,i+107] |> pull()
  
  mod_test[[i]] <- phylo.modularity(mat,partition.gp = labels,
                                    phy=tree_nw,print.progress=T,
                                    CI = TRUE)
}



mod_test |> save(file = 'results/modularity_test.RData')

load('results/modularity_test.RData')
