library(geomorph)
library(tidyverse)
library(broom)
library(OUwie)
library(phytools)
library(ggx)


mat <- read.csv("data/1.processed/Broc2018_maindataset.csv")
clust <- read.csv("data/4.InteractionScores_tSNE/tSNE_Clustermembership_ppx706.csv")

mat1 <- full_join(mat,clust,by="drugdrug")
mat2 <- mat1 %>% 
  filter(!is.na(ebw)&!is.na(ecr)&!is.na(seo)&!is.na(stm)&!is.na(pae)&!is.na(pau))

# Read phylogenetic tree
species <- c("ebw","ecr","stm","seo","pae","pau")
treetype <- "bac_species.txt.hcp"
treefile <- file.path('data/2.Phylogenetics_PhySpeTree',treetype,"iqtree.tree.treefile")
tree_nw <- read.tree(treefile)
# Drop the species that are not needed.
not.species <- tree_nw$tip.label[!(tree_nw$tip.label %in% species)]
tree_nw <- drop.tip(tree_nw,not.species)


# Run test for tSNE clusters
dat <- mat2 %>% select(ebw,ecr,seo,stm,pae,pau) %>% t() 

# Run test for drug categories
RT_tSNE <- compare.multi.evol.rates(A=dat,gp=mat2$clusters,phy=tree_nw,iter=999)
RT_cc <- compare.multi.evol.rates(A=dat,gp=mat2$categorycategory,phy=tree_nw,iter=999)
RT_uu <- compare.multi.evol.rates(A=dat,gp=mat2$useuse,phy=tree_nw,iter=999)
RT_pp <- compare.multi.evol.rates(A=dat,gp=mat2$processprocess,phy=tree_nw,iter=999)
RT_c <- compare.multi.evol.rates(A=dat,gp=mat2$drug_category,phy=tree_nw,iter=999)
RT_p <- compare.multi.evol.rates(A=dat,gp=mat2$targeted_process,phy=tree_nw,iter=999)
RT_u <- compare.multi.evol.rates(A=dat,gp=mat2$use,phy=tree_nw,iter=999)
str(RT_tSNE)
## With network distance.
distance <- read.csv("data/3.Targets_NetworkDistance/DrugTargets4_distance.csv")

mat1 <- full_join(mat,distance,by="drugdrug")
mat2 <- mat1 %>% 
  filter(!is.na(ebw)&!is.na(ecr)&!is.na(seo)&!is.na(stm)&!is.na(pae)&!is.na(pau)) %>% 
  filter(!is.na(Min_distance))

x <- mat2 %>% select(Min_distance,ebw,ecr,stm,seo,pae,pau)
y <- x %>% select(ebw,ecr,stm,seo,pae,pau) %>% t()
RT_dist <- compare.multi.evol.rates(A=y,gp=x$Min_distance,phy=tree_nw,iter=999)

####
RT <- RT_tSNE
df <- data.frame(rates=RT$sigma.d.gp,
           groups=RT$groups)
ggplot(df,aes(x=as.numeric(df$groups),y=df$rates))+
  geom_point()+
  theme_minimal()+
  ggtitle(label = RT$call)+
  theme(axis.text.x = element_text(angle = 90)) 


summary(RT_uu)
plot(RT_tSNE)





# Rates from tSNE clusters ~ distance
rate_tSNE_df <- data.frame(rates=RT_tSNE$sigma.d.gp,
                           clusters=as.numeric(RT_tSNE$groups))

mat_A <- full_join(mat2,clust,by="drugdrug")
mat_B <- full_join(mat_A,rate_tSNE_df,by="clusters")

rat_mat <- mat_B %>% select(Min_distance,rates,clusters) %>% 
  filter(!is.na(Min_distance))

ggplot(rat_mat,aes(x=as.factor(Min_distance),y=rates))+
  geom_point()+
  geom_boxplot()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90)) 

rat_mat %>% group_by(Min_distance) %>% 
  summarise(mean_rate=mean(rates)) %>% 
  ggplot(aes(x=Min_distance,y=mean_rate))+
  geom_point()+
  theme_minimal()
