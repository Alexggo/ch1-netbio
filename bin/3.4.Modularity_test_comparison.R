#!/usr/bin/env Rscript

#' ---
#' title: "3.4.Modularity test comparison"
#' author: "Alex Gil"
#' date: "1/26/2021"
#' output: html_document
#' 
#' 
#' 
## -------------------------------------------------------------------------
library(pacman)
p_load(geomorph,tidyverse,broom,ape,ggx,geomorph,furrr,tictoc,patchwork)
library(ggrepel)


full_dataset <- read_csv("data/3.InteractionScores_tSNE/all_ppx_large_set2.csv")
load('results/set/modularity_tests-115-455_sen.RData') #129 ppx
load(file = 'results/set/rate_tests-5_most_modular.RData')


ptsne.ppx <- seq(5,2505,10) #251 ppx. 
#ALL DDI 175-1455
# SET2 115-455
ptsne.ppx_tested <- ptsne.ppx[ptsne.ppx>=115&ptsne.ppx<=455] #129
names(mod_test) <- ptsne.ppx_tested


labels_df <- full_dataset |> select(paste0("ppx_",ptsne.ppx_tested))
labels <- list()
for (i in 1:dim(labels_df)[2]){
  labels[[i]] <- labels_df[,i] |> pull()
}

number_of_clusters <- labels |> 
  map(unique) |> 
  map_dbl(length)

CR <- c()
CInterval <- list()
P.value <- c()
Z <- c()
CR.mat <- list()


for (i in 1:35){
  CR[i] <- mod_test[[i]]$CR
  CInterval[[i]] <- mod_test[[i]]$CInterval
  P.value[i] <- mod_test[[i]]$P.value
  Z[i] <- mod_test[[i]]$Z
  CR.mat[[i]] <- mod_test[[i]]$CR.mat
}

df <- data.frame(Z_effect=Z,
                 ptsne=ptsne.ppx_tested,
                 p.value=P.value,
                 CR=CR,
                 cluster_n=number_of_clusters) 

df |> 
  filter(p.value<0.05) |> 
  ggplot(aes(x=ptsne,y=CR))+
  geom_point()+
  theme_minimal()+
  xlab("tsne perplexity")

df |> 
  filter(p.value<0.05) |> 
  ggplot(aes(x=ptsne,y=Z_effect))+
  geom_point()+
  theme_minimal()+
  xlab("tsne perplexity")+
  geom_smooth(method="lm")

df |> 
  filter(ptsne<=500) |> 
  ggplot(aes(x=CR,y=Z_effect,label = ptsne))+
  theme_minimal()+
  geom_point()+
  geom_text_repel()+
  geom_hline(yintercept=-26.5)+
  geom_vline(xintercept=0.92)+
  geom_smooth(method="lm")

#SET. 135 ppx

library(broom)
test <- lm(Z_effect~ptsne,df)
broom::augment(test) |> 
  arrange(.resid)

# Adams,2017. The 5 most negative effect size (ZCR) clustering solutions were identified, 
# these represented the hypotheses with the strongest modular signals.

# filtered_df <- df |> 
#   filter(p.value<0.05) |> 
#   filter(Z_effect<=-26.5) |> 
#   filter(CR<0.925) |> 
#   arrange(Z_effect)


filtered_df <- df |> 
  filter(p.value<0.05) |> 
  filter(Z_effect<=-9) |> 
  arrange(Z_effect)

dataset <- full_dataset |> 
  select(1:86,paste0("ppx_",filtered_df$ptsne))

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
# Test differences in rates
dat <- dataset %>%  
  select(ebw,ecr,seo,stm,pae,pau) %>% t() 

# To run test, run the following:
RT_list <- list()
for (i in 1:dim(filtered_df)[1]){
  t_ppx <- filtered_df$ptsne[i]
  t_ppx <- paste0("ppx_",t_ppx)
  groups_t <- dataset |> select(starts_with(t_ppx)) |> pull()
  RT_list[[i]] <- compare.multi.evol.rates(A=dat,gp=groups_t,phy=tree_nw,iter=999)
}
names(RT_list) <- paste0("ppx_",filtered_df$ptsne)
RT_list |> save(file = 'results/set/rate_tests-2_most_modular.RData')

sigma.d.ratio<- c()
p.value_test<- c()
Z_test<- c()
sigma.d.gp <- list()
for (i in 1:length(RT_list)){
  sigma.d.ratio[i] <- RT_list[[i]]$sigma.d.ratio
  p.value_test[i] <- RT_list[[i]]$P.value
  Z_test[i] <- RT_list[[i]]$Z
  sigma.d.gp[[i]] <- RT_list[[i]]$sigma.d.gp
}

df_rate_test <- data.frame(ptsne=filtered_df$ptsne,
                           sigma.d.ratio=sigma.d.ratio,
                           Z_test=Z_test,
                           p.value.rate_test=p.value_test) |> 
  slice_min(p.value.rate_test,n=1)

summary_test <- inner_join(filtered_df,df_rate_test,by="ptsne")

best <- RT_list[["ppx_135"]]  #ALL DDI 245. SET2 PPX 135
summary(best)
plot(best)

RT_result <- data.frame(clusters=as.integer(best$groups),
                        sigma.rate=best$sigma.d.gp)

RT_result %>% ggplot(aes(x=as.factor(clusters),y=sigma.rate))+
  geom_bar(stat="identity")+
  ggtitle("Evolutionary rates for each cluster")+
  xlab("Cluster")+
  ylab(expression("Sigma rate ("~Bliss^2/MYA~")"))+
  theme_minimal()

RT_result$clusters <- factor(RT_result$clusters, levels = RT_result$clusters[order(RT_result$sigma.rate)])
p1 <- ggplot(RT_result, aes(x = clusters, y = sigma.rate)) + 
  geom_bar(stat = "identity")+
  ggtitle("Evolutionary rates for each cluster")+
  xlab("Cluster")+
  ylab(expression("Sigma rate ("~Bliss^2/MYA~")"))+
  theme_minimal()
  


RT_result <- RT_result %>% mutate(clusters=as.numeric(as.character(clusters)))

full_dataset$clusters <- full_dataset$ppx_135 
full_df <- left_join(full_dataset,RT_result,by="clusters")
full_df <- full_df |> select(!starts_with("ppx"))


r <- c()
for (i in 1:dim(full_df)[1]){
  x <- full_df[i,25:30] |> t() |> as.vector()
  x <- x |> unique() |> sort()
  r[i] <- paste0(x,collapse = "-")
}
full_df$type <- r

full_df |> 
  write_csv("data/4.PhylogeneticComparativeMethods/DDI_table_rates_set.csv")
