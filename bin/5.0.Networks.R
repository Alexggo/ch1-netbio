#' ---
#' title: "5.0.Standardize networks"
#' author: "Alex Gil"
#' 
library(igraph)
library(tidyverse)
library(BioNet)
library(ape)
library(geomorph)

# Fix metabolic network.
edges <- read_tsv("data/5.Targets_NetworkDistance/MODEL1108160000_BASE.txt")
node_info <- read_csv("data/5.Targets_NetworkDistance/MODEL1108160000_BASE_node.csv") %>% 
  select(name,GENE_ASSOCIATION) %>% 
  mutate(GENE=gsub("\\(","",as.character(GENE_ASSOCIATION))) %>% 
  mutate(GENE=gsub("\\)","",as.character(GENE))) %>% 
  mutate(GENE=gsub("and","or",as.character(GENE))) %>% 
  separate(GENE,sep="or",into=c("c1","c2","c3","c4",
                                "c5","c6","c7","c8",
                                "c9","c10","c11","c12",
                                "c13","c14","c15","c16",
                                "c17","c18")) %>% 
  select(-GENE_ASSOCIATION) %>% 
  pivot_longer(names_to = "col",values_to = "gene",2:19) %>% 
  drop_na() %>% 
  select(name,gene) %>% 
  mutate(gene=gsub(" ","",as.character(gene))) %>% 
  distinct() %>% 
  rename(node1=name,gene1=gene)

d1 <- left_join(edges,node_info,by="node1")
node_info <- node_info %>% rename(node2=node1,
                                  gene2=gene1)
d2 <- left_join(d1,node_info,by="node2")
d3 <- d2 %>% 
  mutate(gene1=ifelse(is.na(gene1),node1,gene1)) %>% 
  mutate(gene2=ifelse(is.na(gene2),node2,gene2)) %>% 
  select(gene1,gene2) %>% 
  rename(node1=gene1,node2=gene2)
write_tsv(d3,"data/5.Targets_NetworkDistance/MODEL1108160000_edgelist.txt")

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
  #- MODEL1108160000_edgelist.txt

filename <- c("EcoCyc.goldstandardset.txt","EcoliNet.v1.txt","GN.INT.EcoliNet.3568gene.23439link.txt",
              "GO-BP.goldstandardset.txt","CC.EcoliNet.v1.2296gene.50528link.txt","CX.INT.EcoliNet.v1.4039gene.67494link.txt",
              "DC.EcoliNet.2283gene.9643link.txt","EcoliNet.v1.benchmark.txt","HT.INT.EcoliNet.3209gene.15543link.txt",
              "LC.INT.EcoliNet.764gene.1073link.txt","PG.INT.EcoliNet.v1.1817gene.17504link.txt"
              ,"MODEL1108160000_edgelist.txt"
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
  listx[[9]],listx[[10]],listx[[11]],listx[[12]]) %>% unique() %>% sort()

writeLines(keggid_1,file.path("data/5.Targets_NetworkDistance","kegg.ids.txt"))

# Retrieve Uniprot IDs from website
x2 <- read.delim(file.path("data/5.Targets_NetworkDistance","kegg.ids.uniprot.tab"))
colnames(x2)[1] <- "name"
x2 <- x2 %>% select(Entry,name)
x3.node1 <- x2 %>% 
  rename(node1 = name)
x3.node2 <- x2 %>% 
  rename(node2 = name)

#Substitute KEGG_ID codes for Uniprot_ID codes in each network.
for (i in 1:length(filename)){
  t <- file.path("data/5.Targets_NetworkDistance",filename[i])
  x1 <- read_tsv(t)
  x1 <- x1 %>% select(1:2)
  x1$Uni_N1 <- paste0("eco:",x1$node1)
  x1$Uni_N2 <- paste0("eco:",x1$node2)
  x1 <- x1 %>% select(Uni_N1,Uni_N2)
  write.csv(x1,file=file.path("data/5.Targets_NetworkDistance",
                              paste0(filename[i],"net.txt")),row.names = FALSE)
}




