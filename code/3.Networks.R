#' ---
#' title: "5. Network Analysis"
#' author: "Alex Gil"
#' 
library(igraph)
library(tidyverse)

# E.coli interactome
network1 <- read.csv("data/3.Targets_NetworkDistance/2013.1-Ecoli_interactome_clean.csv") %>% 
  select(c(3,4)) %>% as.matrix()
g=graph.edgelist(network1,directed=FALSE)
average.path.length(g)
dist_mat <- shortest.paths(g) #Gives shortest path, more than one may exists between two vertices.

# Known targets
nodes <- read.csv("data/3.Targets_NetworkDistance/DrugTargets3_ecoli.csv")
nodes1 <- nodes %>% select(12) %>% unique()
nodes1 <- nodes1$UNIPROT_eco

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
  filter(Length!=0) %>% 
  filter(Length!=Inf)

# Map Targets to Drugs.
unique_targets <- tabnew2$Uniprot1 %>% unique()
drug_target <-nodes %>% filter(UNIPROT_eco %in% unique_targets)
drug_target1 <- drug_target %>% 
  select(Drug,UNIPROT_eco)
colnames(drug_target1) <- c("Drug1","Uniprot1")

tab1 <- inner_join(tabnew2,drug_target1,by="Uniprot1")
colnames(drug_target1) <- c("Drug2","Uniprot2")
tab1 <- inner_join(tab1,drug_target1,by="Uniprot2")

tab1 <- tab1 %>% filter(tab1$Drug1!=tab1$Drug2)%>% distinct()

y <- tab1 %>% group_by(Drug1,Drug2) %>% 
  summarise(min.length=min(Length))

# Convert to square matrix, to obtain upper diagonal values.
z <- y %>% pivot_wider(names_from=Drug2,values_from=min.length) %>% 
  as.data.frame()
rownames(z) <- z$Drug1
z <- z %>% select(-Drug1)

drug_mat <- z[ , order(names(z))]
drug_mat$D1 <- rownames(drug_mat)
drug_mat <- drug_mat %>%
  arrange(D1) %>% 
  select(-D1)
drug_mat[lower.tri(drug_mat)] <- NA
drug_mat$D1 <- colnames(drug_mat)
drug_mat <- drug_mat %>%
  select(dim(drug_mat)[2],1:(dim(drug_mat)[2]-1))

drug_mat1 <- drug_mat  %>% 
  pivot_longer(names_to="D2",values_to="Min_distance",2:(dim(drug_mat)[2]))%>%
  filter(!is.na(Min_distance)) %>% 
  filter(Min_distance!=0) 
drug_mat2 <- drug_mat1 %>% mutate(drugdrug=paste0(D1,"-",D2))

write.csv(drug_mat2,"data/3.Targets_NetworkDistance/DrugTargets4_distance.csv",row.names = F)


### Compare interaction scores and rates to network distance between targets.
# Interaction score
scores <- read.csv("data/1.processed/Broc2018_maindataset.csv")
x1 <- scores %>% select(drug_pair,ebw,ecr,BM.sigmasq,OU.alpha)

graph1 <- full_join(drug_mat2,x1,by="drug_pair")

ggplot(graph1) +
  aes(x = Mean_length, y = OU.alpha) +
  geom_point(size = 1L, colour = "#0c4c8a") +
  theme_minimal()+
  geom_smooth(method = "lm", se = FALSE)
  
