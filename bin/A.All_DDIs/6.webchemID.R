#MAPPING

library(igraph)
library(tidyverse)
library(BioNet)
library(ape)
library(geomorph)

filename <- c("EcoCyc.goldstandardset.txt","EcoliNet.v1.txt","GN.INT.EcoliNet.3568gene.23439link.txt",
              "GO-BP.goldstandardset.txt","CC.EcoliNet.v1.2296gene.50528link.txt","CX.INT.EcoliNet.v1.4039gene.67494link.txt",
              "DC.EcoliNet.2283gene.9643link.txt","EcoliNet.v1.benchmark.txt","HT.INT.EcoliNet.3209gene.15543link.txt",
              "LC.INT.EcoliNet.764gene.1073link.txt","PG.INT.EcoliNet.v1.1817gene.17504link.txt"
              ,"MODEL1108160000_edgelist.txt")


# Map drugs to targets (for those that are known)
nodes <- read.csv("data/5.Targets_NetworkDistance/DrugTargets3_ecoli.csv")
nodes1 <- nodes %>%
  select(KEGG_eco)%>%
  unique() %>%
  pull()

t <- file.path("data/5.Targets_NetworkDistance",paste0(filename,"net.txt"))
x1 <- map(t,read_csv)
x1 <- x1 %>% map(select,Uni_N1,Uni_N2)
g1 <- x1 %>% map(graph_from_data_frame, directed = FALSE, vertices = NULL)
deg <- g1 %>% map(igraph::degree)

degree_distribution(g1) %>% plot()

g1 %>% map(map(igraph::degree_distribution),.) %>% plot()
