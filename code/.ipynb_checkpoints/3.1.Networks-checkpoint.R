# Do proteins with annotations share the same rate than proteins without annotations.
#1. Use network to get groups
#2. Look at rate in tSNE
#3. Are they in the same tSNE group.

library(igraph)
library(tidyverse)
library(BioNet)
library(ape)
library(geomorph)


# Load data and tree
mat <- read.csv("data/1.processed/Broc2018_maindataset.csv")
clust <- read.csv("data/4.InteractionScores_tSNE/tSNE_Clustermembership_ppx706.csv")
mat1 <- full_join(mat,clust,by="drugdrug")
mat2 <- mat1 %>% 
  filter(!is.na(ebw)&!is.na(ecr)&!is.na(seo)&!is.na(stm)&!is.na(pae)&!is.na(pau))
species <- c("ebw","ecr","stm","seo","pae","pau")
treetype <- "bac_species.txt.hcp"
treefile <- file.path('data/2.Phylogenetics_PhySpeTree',treetype,"iqtree.tree.treefile")
tree_nw <- read.tree(treefile)
not.species <- tree_nw$tip.label[!(tree_nw$tip.label %in% species)]
tree_nw <- drop.tip(tree_nw,not.species)


# Load network. 
#- EcoCyc.goldstandardset.txt
#- EcoliNet.v1.benchmark.txt
#- HT.INT.EcoliNet.3209gene.15543link.txt
#- LC.INT.EcoliNet.764gene.1073link.txt

filename <- "HT.INT.EcoliNet.3209gene.15543link.txt"
x6 <- read.csv(file.path("data/3.Targets_NetworkDistance",paste0(filename,"network.tab")))
g1 <- graph_from_data_frame(x6, directed = FALSE, vertices = NULL)

# Relate target interactions to drug drug interactions.

# Known targets
nodes <- read.csv("data/3.Targets_NetworkDistance/DrugTargets3_ecoli.csv")
nodes1 <- nodes %>%
  select(UNIPROT_eco)%>%
  unique() %>%
  pull()

adj <- as_adjacency_matrix(g1, type = "both",names = TRUE,
                           sparse = igraph_opt("sparsematrices")) %>% as.matrix()

# Filter columns and rows for known targets
r <- rownames(adj) %in% nodes1
c <- colnames(adj) %in% nodes1
filtered <- adj[r,c] %>% 
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
  pivot_longer(names_to="Uniprot2",values_to="value",2:(dim(tabnew)[2]))%>%
  filter(!is.na(value)) 

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

y <- tab1 %>% group_by(Drug1,Drug2)  %>% 
  distinct() 
y$DDI <- NA

for (i in 1:dim(y)[1]){
  nam <- sort(c(y$Drug1[i],y$Drug2[i]),decreasing = F)
  
  y$DDI[i] <- paste0(nam, collapse="_")
}

y1 <- y %>% select(DDI,value) %>% distinct() %>% 
  mutate(Drug_pair=DDI)%>% 
  separate(col=DDI,into=c("Drug1","Drug2"),sep="_")%>% 
  mutate(connection=ifelse(value==TRUE,"yes","no"))


# Merge both datasets.
mat3 <- mat2 %>% mutate(Drug_pair=drug_pair.x)
merged <- inner_join(mat3,y1,by="Drug_pair") %>% 
  select(Drug_pair,Drug1, Drug2,connection,clusters)

# Contingency tables and chi squared test.
# Problem: every DDI belong to only one cluster. So I can not test same/different cluster vs yes/no connection.
# (Ask Josh)
# Do DDI that are part of the same cluster have connections in the network? 
ggplot(merged,aes(x=clusters,fill=connection))+
  geom_bar(position = "fill")+
  theme_minimal()

# For DDI with targets with connections (yes/no)
# Compare tSNE rate distributions barchart.
dat <- mat2 %>% select(ebw,ecr,seo,stm,pae,pau) %>% t() 
RT_tSNE <- compare.multi.evol.rates(A=dat,gp=mat2$clusters,phy=tree_nw,iter=999)

RT <- RT_tSNE
df <- data.frame(rates=RT$sigma.d.gp,
                 groups=RT$groups)
ggplot(df,aes(x=as.numeric(df$groups),y=df$rates))+
  geom_point()+
  theme_minimal()+
  ggtitle(label = RT$call)+
  theme(axis.text.x = element_text(angle = 90)) 

summary(RT_tSNE)
plot(RT_tSNE)


df <- df %>% mutate(clusters=as.numeric(groups)) %>% select(-groups)
new_df <- inner_join(merged,df,by="clusters")



# Are differences in means significant?
yes1 <- new_df %>% filter(connection=="yes") %>% select(rates) %>% pull()
no1 <- new_df %>% filter(connection=="no") %>% select(rates) %>% pull()
t <- t.test(x=no1,y=yes1) #Two-sided

ggplot(new_df,aes(x=connection,y=rates,fill=connection))+
  geom_boxplot()+
  geom_jitter()+
  theme_minimal()


# Classify target relations by connectivity level.
# bins: low/medium/high connectivity (x axis).
# y-axis: 1. average rate from tSNE cluster 


conn <- read.csv(file.path("data/3.Targets_NetworkDistance",paste0(filename,"connectmat.csv")))
hist(conn$value)


conn <- conn %>% mutate(Drug_pair=paste0(Drug1,"_",Drug2))

new_df <- inner_join(mat3,conn,by="Drug_pair") %>% 
  select(Drug_pair,Drug1,Drug2,value,clusters)

dat1 <- inner_join(new_df,df,by="clusters")

ggplot(dat1,aes(x=value,y=rates))+
  geom_point()+
  theme_minimal()+
  xlab("Connectivity")

# bins: low/medium/high connectivity (x axis).
# 2. % that are in the same cluster. This is wrong (each DDI is in only one cluster).




# Delta-rate (proteins) ~ PPI(+/-) or distance/connectivity. Ask Josh again for clarification.







