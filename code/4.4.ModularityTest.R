library(geomorph)
library(tidyverse)
library(broom)
library(ape)

mat <- read.csv("data/1.processed/Broc2018_maindataset.csv") %>%
  filter(!is.na(ebw)&!is.na(ecr)&!is.na(seo)&!is.na(stm)&!is.na(pae)&!is.na(pau))

# Read phylogenetic tree
species <- c("ebw","ecr","stm","seo","pae","pau")

treetype <- "bac_species.txt.hcp"
treefile <- file.path('data/2.Phylogenetics_PhySpeTree',treetype,"iqtree.tree.treefile")
tree_nw <- read.tree(treefile)
# Drop the species that are not needed.
not.species <- tree_nw$tip.label[!(tree_nw$tip.label %in% species)]
tree_nw <- drop.tip(tree_nw,not.species)

# 13 Clusters
clusters <- read.csv("data/4.InteractionScores_tSNE/tSNE_Clustermembership_ppx706.csv")

mat1 <- full_join(mat,clusters,by="drugdrug")
mat2 <- mat1[,12:17] %>% t()


mod1 <- phylo.modularity(mat2,partition.gp = mat1$clusters,phy=tree_nw,print.progress=T,
                         CI = TRUE)
mod1
# CR: 0.9231
# P-value: 0.001
# Effect Size: -26.3594
# Based on 1000 random permutations
# Confidence Intervals 0.8643 
# Confidence Intervals 1.0121


# 4 Clusters
clusters <- read.csv("data/4.InteractionScores_tSNE/tSNE_Clustermembership_ppx825.csv")

mat1 <- full_join(mat,clusters,by="drugdrug")
mat2 <- mat1[,12:17] %>% t()


mod2 <- phylo.modularity(mat2,partition.gp = mat1$clusters,phy=tree_nw,print.progress=T,
                         CI = TRUE)

mod2
# CR: 0.985
# P-value: 0.001
# Effect Size: -13.9241
# Based on 1000 random permutations
#Confidence Intervals 0.9544 
#Confidence Intervals 1.0049


# Adams,2017. The most negative effect size (ZCR) is identified, 
# this represented the hypothesis representing the strongest modular signal.

#In this case model 1 is better.

