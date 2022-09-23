library(pacman)
p_load(broom,tidyverse,OUwie,phytools,ape,ggx,geomorph)

mat <- read.csv("data/1.processed/Broc2018_maindataset.csv")
clust <- read.csv("data/3.InteractionScores_tSNE/tSNE_Clustermembership_ppx706.csv")

mat1 <- full_join(mat,clust,by="drug_pair")
mat2 <- mat1 %>% 
  filter(!is.na(ebw)&!is.na(ecr)&!is.na(seo)&!is.na(stm)&!is.na(pae)&!is.na(pau))

# Read phylogenetic tree
species <- c("Escherichia_coli_K-12_ebw",
             "Escherichia_coli_O8_IAI1_ecr",
             "Salmonella_enterica_serovar_Typhimurium_LT2_stm",
             "Salmonella_enterica_serovar_Typhimurium_14028S_seo",
             "Pseudomonas_aeruginosa_PAO1_pae",
             "Pseudomonas_aeruginosa_UCBPP-PA14_pau")

treetype <- "bac_species.txt.hcp"
treefile <- file.path('data/2.2.Phylogenetics_Bayesian/DDI_BD_str',"concatenate.trees.tre")
tree_nw <- read.nexus(treefile)
# Drop the species that are not needed.
not.species <- tree_nw$tip.label[!(tree_nw$tip.label %in% species)]
tree_nw <- drop.tip(tree_nw,not.species)
tree_nw$tip.label <- c("ebw","ecr","pae","pau","seo","stm")


# Run test for drug categories
dat <- mat2 %>%  
  select(ebw,ecr,seo,stm,pae,pau) %>% t() 
RT_tSNE <- compare.multi.evol.rates(A=dat,gp=mat2$clusters,phy=tree_nw,iter=999)

#RT_tSNE is significant.
summary(RT_tSNE)
plot(RT_tSNE)

RT_result <- data.frame(clusters=as.integer(RT_tSNE$groups),
                        sigma.rate=RT_tSNE$sigma.d.gp)

RT_result %>% ggplot(aes(x=as.factor(clusters),y=sigma.rate))+
  geom_bar(stat="identity")+
  ggtitle("Evolutionary rates for each cluster")+
  xlab("Cluster")+
  ylab("Sigma rate")

RT_result$clusters <- factor(RT_result$clusters, levels = RT_result$clusters[order(RT_result$sigma.rate)])
ggplot(RT_result, aes(x = clusters, y = sigma.rate)) + theme_bw() + 
  geom_bar(stat = "identity")+
  xlab("Cluster")+
  ylab("Sigma rate")+
  ggtitle("Evolutionary rates for each cluster")+
  xlab("Cluster")+
  ylab("Sigma rate")

RT_result <- RT_result %>% mutate(clusters=as.numeric(as.character(clusters)))
x <- inner_join(mat2,RT_result,by="clusters") |> 
  mutate(Drug1=Drug1.x) %>% 
  mutate(Drug2=Drug2.x) %>% 
  select(-c(Drug1.x,Drug2.x,Drug1.y,Drug2.y))

write.csv(x,"data/4.PhylogeneticComparativeMethods/dataset1_ratios.csv",row.names = FALSE)

#RT_tSNE is not significant if cluster 13 is excluded.
dat <- mat2 %>% filter(clusters!=13) %>%
 select(ebw,ecr,seo,stm,pae,pau) %>% t()
tags <- mat2 %>% filter(clusters!=13) %>%
 select(clusters) %>% pull()
RT_tSNE <- compare.multi.evol.rates(A=dat,gp=tags,phy=tree_nw,iter=999)
summary(RT_tSNE)
plot(RT_tSNE)
