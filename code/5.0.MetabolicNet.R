#' ---
#' title: "5.0.MetabolicNet"
#' author: "Alex Gil"
#' aim: Convert metabolic stoichiometric matrix into network
library(igraph)
library(tidyverse)
library(BioNet)


filename <- "ecoli_core_model_Smatrix.csv"
t <- file.path("data/3.Targets_NetworkDistance",filename)
x1 <- read.csv(t,sep=",") 
rownames(x1) <- x1$ï..METABOLITE...ABBREVIATION
x1 <- x1 %>% select(-1) %>% t()

# X1 is the incidence matrix of a bipartite graph, we want to transform it into an adjacency matrix keeping only the enzymes.
g1 <- graph_from_incidence_matrix(x1, directed = FALSE,
                            mode = "total")

mat1 <- as_edgelist(g1,names=TRUE)

filename <- "ecoli_core_model_reactions.csv"
t <- file.path("data/5.Targets_NetworkDistance",filename)
x2 <- read.csv(t,sep=",") 
colnames(x2)[1]<- "Name"

x3 <- x2 %>% select(Name,geneAssociation) %>% 
  filter(!geneAssociation=="")



colnames(mat1)<- c("N1","N2")
d1 <- mat1 %>% as.data.frame()
x4 <- x3
colnames(x4) <- c("N1","G1")

a <- full_join(d1,x4,by="N1")

b <- separate(a,col="G1",into=paste0("A",1:80),sep=" ")

tab <- vector()
for (i in 1:dim(b)[1]){
  vector <- b[i,3:82] 
  logic <- !(vector %in% c(""," ","(",")","or","and"))
  tab[i] <- vector[logic] %>% sort() %>% paste0(collapse = "-")
  
}
b$A <- tab
c <- b %>% select(N1,N2,A) 
d <- separate(c,col="A",into=paste0("A",1:17),sep="-")
e <- d %>% pivot_longer(names_to = "colname",values_to = "geneID",3:19) %>% 
  select(-colname) %>% distinct() %>% filter(!is.na(geneID)) %>% 
  filter(!is.na(N2)) %>% filter(N1!="Biomass_Ecoli_core")

colnames(e) <- c("genename_N1","node2","node1")
e <- e %>% select(node1,node2,genename_N1)


write.table(e,file.path("data/5.Targets_NetworkDistance",paste0(filename,"net.txt")),row.names = FALSE,sep="\t")
