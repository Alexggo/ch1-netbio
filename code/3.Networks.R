#' ---
#' title: "5. Network Analysis"
#' author: "Alex Gil"
#' 
library(igraph)
library(tidyverse)
library(BioNet)

#First fix the PPI network from E.coli adding KEGG ids.
  #- EcoCyc.goldstandardset.txt
  #- EcoliNet.v1.benchmark.txt
  #- HT.INT.EcoliNet.3209gene.15543link.txt
  #- LC.INT.EcoliNet.764gene.1073link.txt

filename <- "LC.INT.EcoliNet.764gene.1073link.txt"
t <- file.path("data/3.Targets_NetworkDistance",filename)
x1 <- read_tsv(t) 
x1$node1 <- paste0("eco:",x1$node1)
x1$node2 <- paste0("eco:",x1$node2)
keggid_1 <- c(x1$node1,x1$node2) %>% unique()
writeLines(keggid_1,file.path("data/3.Targets_NetworkDistance",paste0(filename,".link_ids.txt")))
           
#Retrieve Uniprot IDs.
x2 <- read.delim(file.path("data/3.Targets_NetworkDistance",paste0(filename,".link_ids.txt","_uniprot.tab")))
                           
colnames(x2)[1] <- "name"
x3 <- x2 %>% mutate(node1=name) %>% select(node1,Entry)

#Substitute KEGG_ID network for Uniprot_ID network.
x4<- inner_join(x1,x3,by="node1")
x4$N1 <- x4$Entry
x4 <- x4 %>% select(node1,N1,node2)

colnames(x3)[1] <- "node2"
x4<- inner_join(x3,x4,by="node2")
x4$N2 <- x4$Entry
x5 <- x4 %>% select(N1,node1,N2,node2)
x6 <- x5 %>% select(N1,N2) %>%  distinct()

write.csv(x6,file.path("data/3.Targets_NetworkDistance",paste0(filename,"network.tab")),row.names = F)

#Graph
g1 <- graph_from_data_frame(x6, directed = FALSE, vertices = NULL)


#Calculate average path length in the network.
average.path.length(g1)
#Calculate the shortest path between nodes.
dist_mat <- shortest.paths(g1) 

# Known targets
nodes <- read.csv("data/3.Targets_NetworkDistance/DrugTargets3_ecoli.csv")
nodes1 <- nodes %>%
  select(UNIPROT_eco)%>%
  unique() %>%
  pull()

# Filter columns and rows for known targets
r <- rownames(dist_mat) %in% nodes1
c <- colnames(dist_mat) %in% nodes1
filtered <- dist_mat[r,c] %>% 
  as.data.frame()
tabnew <- filtered
# Convert matrix to table
tabnew <- tabnew[ , order(names(tabnew))]
tabnew$Uniprot1 <- rownames(tabnew)
tabnew <- tabnew %>%
  arrange(Uniprot1) %>% 
  select(-Uniprot1)

tabnew[lower.tri(tabnew)] <- NA
tabnew$Uniprot1 <- colnames(tabnew)
tabnew1 <- tabnew %>%
  select(dim(tabnew)[2],1:(dim(tabnew)[2]-1))
tabnew2 <- tabnew1  %>% 
  pivot_longer(names_to="Uniprot2",values_to="Length",2:(dim(tabnew)[2]))%>%
  filter(!is.na(Length)) %>% 
  filter(Length!=0)

# Map Targets to Drugs.
unique_targets <- tabnew2$Uniprot1 %>% unique()
drug_target <-nodes %>%
  filter(UNIPROT_eco %in% unique_targets)
drug_target1 <- drug_target %>% 
  select(Drug,UNIPROT_eco)
colnames(drug_target1) <- c("Drug1","Uniprot1")


tab1 <- inner_join(tabnew2,drug_target1,by="Uniprot1")
colnames(drug_target1) <- c("Drug2","Uniprot2")
tab1 <- inner_join(tab1,drug_target1,by="Uniprot2")

tab1 <- tab1 %>% filter(tab1$Drug1!=tab1$Drug2)%>% distinct()

y <- tab1 %>% group_by(Drug1,Drug2) %>% 
  summarise(value=round(mean(Length))) %>% 
  distinct() %>% arrange(value) 
y$DDI <- NA

for (i in 1:dim(y)[1]){
 nam <- sort(c(y$Drug1[i],y$Drug2[i]),decreasing = F)
 
 y$DDI[i] <- paste0(nam, collapse="-")
}

y1 <- y %>% select(DDI,value) %>% distinct() %>% 
  separate(col=DDI,into=c("Drug1","Drug2"),sep="-")

write.csv(y1,file.path("data/3.Targets_NetworkDistance",paste0(filename,"distmat.csv")),row.names = F)

          
          

#Calculate k-edge-connectivity between targets.

#Which nodes are in g1?
filt <- nodes1 %in% rownames(dist_mat)
selection <- nodes1[filt]

t <- matrix(nrow = length(selection),ncol=length(selection))
for (i in 1:length(selection)){
  for (j in 1:length(selection)){
    if(i==j){t[i,j] <- NA}else{
      t[i,j] <- edge_connectivity(g1,source = selection[i],target = selection[j])
    }
  }
}

colnames(t) <- selection
rownames(t) <- selection

# Convert matrix to table
tabnew <- t %>% as.data.frame()
tabnew <- tabnew[ , order(names(tabnew))]
tabnew$Uniprot1 <- rownames(tabnew)
tabnew <- tabnew %>%
  arrange(Uniprot1) %>% 
  select(-Uniprot1)
tabnew[lower.tri(tabnew)] <- NA
tabnew$Uniprot1 <- colnames(tabnew)
tabnew1 <- tabnew %>%
  select(dim(tabnew)[2],1:(dim(tabnew)[2]-1))
tabnew2 <- tabnew1  %>% 
  pivot_longer(names_to="Uniprot2",values_to="K.edge",2:(dim(tabnew)[2]))%>%
  filter(!is.na(K.edge)) %>% 
  filter(K.edge!=0) 

# Map Targets to Drugs.
unique_targets <- tabnew2$Uniprot1 %>% unique()
drug_target <-nodes %>%
  filter(UNIPROT_eco %in% unique_targets)
drug_target1 <- drug_target %>% 
  select(Drug,UNIPROT_eco)
colnames(drug_target1) <- c("Drug1","Uniprot1")


tab1 <- inner_join(tabnew2,drug_target1,by="Uniprot1")
colnames(drug_target1) <- c("Drug2","Uniprot2")
tab1 <- inner_join(tab1,drug_target1,by="Uniprot2")

tab1 <- tab1 %>% filter(tab1$Drug1!=tab1$Drug2)%>% distinct()

y <- tab1 %>% group_by(Drug1,Drug2) %>% 
  summarise(value=round(mean(K.edge))) %>% 
  distinct() %>% arrange(value) 
y$DDI <- NA

for (i in 1:dim(y)[1]){
  nam <- sort(c(y$Drug1[i],y$Drug2[i]),decreasing = F)
  
  y$DDI[i] <- paste0(nam, collapse="-")
}

y1 <- y %>% select(DDI,value) %>% distinct() %>% 
  separate(col=DDI,into=c("Drug1","Drug2"),sep="-")

write.csv(y1,file.path("data/3.Targets_NetworkDistance",paste0(filename,"connectmat.csv")),row.names = F)






