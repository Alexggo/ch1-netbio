# From networks to Uniprot Ids:
library(tidyverse)
eco <- read.table("data/3.Targets_NetworkDistance/STRINGDB/STRINGDB-E.coliK12.v12.0.txt",header = TRUE) |> 
  select(1,2)
ecr <- read.table("data/3.Targets_NetworkDistance/STRINGDB/STRINGDB-EcoliIAI1.v12.0.txt",header = TRUE)|> 
  select(1,2)
pae <- read.table("data/3.Targets_NetworkDistance/STRINGDB/STRINGDB-PAO1.v12.0.txt",header = TRUE)|> 
  select(1,2)
stm <- read.table("data/3.Targets_NetworkDistance/STRINGDB/STRINGDB-STLT2.protein.links.v12.0.txt",header = TRUE)|> 
  select(1,2)

eco_vec <- c(eco$protein1,eco$protein2) |> unique()
ecr_vec <- c(ecr$protein1,ecr$protein2) |> unique()
pae_vec <- c(pae$protein1,pae$protein2) |> unique()
stm_vec <- c(stm$protein1,stm$protein2) |> unique()

kegg_ids <- c(eco_vec,ecr_vec,stm_vec,pae_vec)
write.table(kegg_ids,"data/3.Targets_NetworkDistance/STRINGDB/keggids.txt",
            col.names = FALSE,row.names=FALSE,quote = FALSE)
