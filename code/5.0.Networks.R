#' ---
#' title: "5.0.Standardize networks"
#' author: "Alex Gil"
#' 
library(plyr)
library(igraph)
library(tidyverse)
library(BioNet)
library(ape)
library(geomorph)

#Standardize networks to Uniprot list format.
  #- EcoCyc.goldstandardset.txt
  #- EcoliNet.v1.txt
  #- GN.INT.EcoliNet.3568gene.23439link.txt
  #- GO-BP.goldstandardset.txt
  #- CC.EcoliNet.v1.2296gene.50528link.txt
  #- CX.INT.EcoliNet.v1.4039gene.67494link.txt
  #- DC.EcoliNet.2283gene.9643link.txt
  #- EcoliNet.v1.benchmark.txt
  #- HT.INT.EcoliNet.3209gene.15543link.txt
  #- LC.INT.EcoliNet.764gene.1073link.txt
  #- PG.INT.EcoliNet.v1.1817gene.17504link.txt
  #- ecoli_core_model_reactions.csvnet.txt

filename <- c("EcoCyc.goldstandardset.txt","EcoliNet.v1.txt","GN.INT.EcoliNet.3568gene.23439link.txt",
              "GO-BP.goldstandardset.txt","CC.EcoliNet.v1.2296gene.50528link.txt","CX.INT.EcoliNet.v1.4039gene.67494link.txt",
              "DC.EcoliNet.2283gene.9643link.txt","EcoliNet.v1.benchmark.txt","HT.INT.EcoliNet.3209gene.15543link.txt",
              "LC.INT.EcoliNet.764gene.1073link.txt","PG.INT.EcoliNet.v1.1817gene.17504link.txt"
              #,"ecoli_core_model_reactions.csvnet.txt"
)

listx <- list()
for (i in 1:length(filename)){
  t <- file.path("data/5.Targets_NetworkDistance",filename[i])
  x1 <- read_tsv(t)
  x1$node1 <- paste0("eco:",x1$node1)
  x1$node2 <- paste0("eco:",x1$node2)
  listx[[i]] <- c(x1$node1,x1$node2) %>% unique()
}
keggid_1 <- c(listx[[1]],listx[[2]],listx[[3]],listx[[4]],listx[[5]],listx[[6]],listx[[7]],listx[[8]],
  listx[[9]],listx[[10]],listx[[11]]) %>% unique() %>% sort()

writeLines(keggid_1,file.path("data/5.Targets_NetworkDistance","kegg.ids.txt"))

# Retrieve Uniprot IDs from website
x2 <- read.delim(file.path("data/5.Targets_NetworkDistance","kegg.ids.uniprot.tab"))
colnames(x2)[1] <- "name"
x2 <- x2 %>% select(Entry,name)
x3.node1 <- x2 %>% rename(node1=name)
x3.node2 <- x2 %>% rename(node2=name)

#Substitute KEGG_ID codes for Uniprot_ID codes in each network.
for (i in 1:length(filename)){
  t <- file.path("data/5.Targets_NetworkDistance",filename[i])
  x1 <- read_tsv(t)
  x1 <- x1 %>% select(1:2)
  x1$node1 <- paste0("eco:",x1$node1)
  x1$node2 <- paste0("eco:",x1$node2)
  x1 <- merge(x1, x3.node1, by = 'node1')
  colnames(x1)[3] <- "Uni_N1"
  x1 <- merge(x1, x3.node2, by = 'node2')
  colnames(x1)[4] <- "Uni_N2"
  x1 <- x1 %>% select(node1,node2,Uni_N1,Uni_N2) %>% 
    distinct()
  write.csv(x1,file=file.path("data/5.Targets_NetworkDistance",
                              paste0(filename[i],"net.txt")),row.names = FALSE)
}




