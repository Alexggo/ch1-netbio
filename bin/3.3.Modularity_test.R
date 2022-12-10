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
p_load(geomorph,tidyverse,broom,ape,ggx,geomorph,furrr,tictoc)

availableCores()
plan("multisession",workers=8)
#plan("multicore",workers=8)


full_dataset <- read_csv("data/3.InteractionScores_tSNE/all_ppx_large_set2.csv")

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

# ALL_DDI. Good tsnes between 175-1455. Indexes 18-146
# SET2. Good tsnes between 115-455. Indexes 12-46

mod_test <- list()

labels <- list()
for (i in 1:35){
  labels[[i]] <- full_dataset[,i+97] |> pull()
}

number_of_clusters <- labels |> map(unique) |> map_dbl(length)

test_mod <- function(g){
  mod <- phylo.modularity(mat,partition.gp = g,
                   phy=tree_nw,print.progress=T,
                   CI = TRUE)
  

}

tic()
mod_test <- labels |>
  lapply(test_mod)
toc()

mod_test |> save(file = 'results/set/modularity_tests-115-455_sen.RData')
