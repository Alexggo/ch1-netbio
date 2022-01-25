#' ---
#' title: "5.2.Do connected proteins have faster or slower rates?"
#' author: "Alex Gil"
#' 
library(tidyverse)
library(broom)
library(igraph)
library(patchwork)
library("matrixStats")
library("ComplexHeatmap")
library(ggpubr)
library(plotrix)

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



rates <- "dataset1_ratios.csv"
rates <- read.csv(file.path("data/4.PhylogeneticComparativeMethods",rates))

list1 <- list()
for (i in 1:11){
  # Load each network
  t <- file.path("data/5.Targets_NetworkDistance",paste0(filename[i],"net.txt"))
  x1 <- read_csv(t)
  x1 <- x1 %>% select(Uni_N1,Uni_N2)
  # Get adjacency data from graph
  g1 <- graph_from_data_frame(x1, directed = FALSE, vertices = NULL)
  adj <- as_adjacency_matrix(g1, type = "both",names = TRUE,
                             sparse = igraph_opt("sparsematrices")) %>% 
    as.matrix() %>% 
    as.data.frame()
  distnew2 <- adj[,order(colnames(adj))]
  distnew2$Uniprot1 <- rownames(distnew2)
  distnew2 <- distnew2 %>%
    arrange(Uniprot1) %>% 
    select(-Uniprot1)
  distnew2[lower.tri(distnew2)] <- NA
  distnew2$Uniprot1 <- colnames(distnew2)
  distnew2 <- distnew2 %>%
    select(dim(distnew2)[2],1:(dim(distnew2)[2]-1))
  distnew3 <- distnew2  %>% 
    pivot_longer(names_to="Uniprot2",values_to="adjancency",2:(dim(distnew2)[2]))%>%
    filter(!is.na(adjancency)) %>% 
    mutate(U_ID=paste0(Uniprot1,"-",Uniprot2))
  # Drug targets.
  nodes <- read.csv(file.path("data/5.Targets_NetworkDistance",paste0(filename[i],"netvalues_target.csv")))
  nodes <- nodes %>% 
    mutate(U_ID=paste0(Uniprot1,"-",Uniprot2))
  
  a <- inner_join(distnew3,nodes,by="U_ID")
  # Add rates.
  r <- rates %>% 
    select(clusters,sigma.rate) %>%
    distinct() %>% 
    arrange(sigma.rate)
  r$clusters <- factor(r$clusters, levels = r$clusters[order(r$sigma.rate)])
  
  tot <- inner_join(a,rates,by="drug_pair") %>% 
    distinct() %>% mutate(clusters=as.factor(clusters))
  
  tot$clusters <- factor(tot$clusters, levels = r$clusters)
  
  list1[[i]] <- tot
}
# List 1 contains the complete dataset for each network.
# df_tot contains all networks combined.
df_tot <- bind_rows(list1, .id = "id") %>% 
  mutate(adjancency=as.factor(adjancency)) %>% 
  mutate(connection_groups=cut_interval(K.edge, 6)) %>% 
  mutate(connection_groups=ifelse(is.na(connection_groups),0,connection_groups)) %>% 
  mutate(connection_onoff=ifelse(is.na(K.edge),"disconnected","connected"))  %>% 
  mutate(network=ifelse(network=="EcoCyc.goldstandardset.txt","Co-functional (EcoCyc)",
                         ifelse(network=="EcoliNet.v1.txt","EN: all networks",
                                ifelse(network=="GN.INT.EcoliNet.3568gene.23439link.txt","Similar genomic context",
                                       ifelse(network=="GO-BP.goldstandardset.txt","Co-functional (GO-BP)",
                                              ifelse(network=="CC.EcoliNet.v1.2296gene.50528link.txt","Co-citation",
                                                     ifelse(network=="CX.INT.EcoliNet.v1.4039gene.67494link.txt","Co-expression",
                                                            ifelse(network=="DC.EcoliNet.2283gene.9643link.txt","Co-occurence of prot. domains",
                                                                   ifelse(network=="EcoliNet.v1.benchmark.txt","Co-functional (EcoCyc/GO-BP)",
                                                                          ifelse(network=="HT.INT.EcoliNet.3209gene.15543link.txt","High-throughput PPI",
                                                                                 ifelse(network=="LC.INT.EcoliNet.764gene.1073link.txt","Small/medium-scale PPI",
                                                                                        ifelse(network=="PG.INT.EcoliNet.v1.1817gene.17504link.txt","Similar phylogenetic profiles","NA")))))))))))) %>% 
  mutate(int_sign_ebw = factor(int_sign_ebw,levels = c("Synergy","Additivity","Antagonism"))) %>% 
  mutate(drug_cat=ifelse(drug_category=="Same",1,0)) %>% 
  mutate(targ=ifelse(targeted_process=="Same",1,0)) %>% 
  mutate(us=ifelse(use=="Same",1,0)) %>% 
  mutate(ebw_g=ifelse(int_sign_ebw=="Antagonism",+1,
                      ifelse(int_sign_ebw=="Additivity",0,
                             ifelse(int_sign_ebw=="Synergy",-1,NA)))) %>% 
  mutate(ecr_g=ifelse(int_sign_ecr=="Antagonism",+1,
                      ifelse(int_sign_ecr=="Additivity",0,
                             ifelse(int_sign_ecr=="Synergy",-1,NA)))) %>% 
  mutate(seo_g=ifelse(int_sign_seo=="Antagonism",+1,
                      ifelse(int_sign_seo=="Additivity",0,
                             ifelse(int_sign_seo=="Synergy",-1,NA)))) %>% 
  mutate(stm_g=ifelse(int_sign_stm=="Antagonism",+1,
                      ifelse(int_sign_stm=="Additivity",0,
                             ifelse(int_sign_stm=="Synergy",-1,NA)))) %>% 
  mutate(pae_g=ifelse(int_sign_pae=="Antagonism",+1,
                      ifelse(int_sign_pae=="Additivity",0,
                             ifelse(int_sign_pae=="Synergy",-1,NA)))) %>% 
  mutate(pau_g=ifelse(int_sign_pau=="Antagonism",+1,
                      ifelse(int_sign_pau=="Additivity",0,
                             ifelse(int_sign_pau=="Synergy",-1,NA)))) %>% 
  rowwise() %>% 
  mutate(sum_g=ebw_g+ecr_g+seo_g+stm_g+pae_g+pau_g)  
  

# Relation between minimum distance and interaction type (ebw). 
my_comparisons <- list( c("Synergy", "Additivity"), c("Additivity", "Antagonism"), c("Synergy", "Antagonism") )
net <- df_tot$network %>% unique()
listA <- list()
for (i in 1:11){
  
  font_size <- 2
  listA[[i]] <- df_tot %>% filter(network==net[i]) %>% 
    ggboxplot(x = "int_sign_ebw", y = "path.length",
              fill = "int_sign_ebw")+ 
    stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = -3,size=font_size) +
    stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
    theme_minimal()+
    scale_fill_brewer(palette="Set2")+
    xlab("Type of interaction")+
    ylab("Minimum path length between proteins")+
    theme(legend.position = "none")+
    ggtitle(net[i])
  
  listA[[i]]$layers[[2]]$aes_params$textsize <- font_size
  
}
wrap_plots(listA[[1]],listA[[2]],listA[[3]],listA[[4]],listA[[5]],listA[[6]],
           listA[[7]],listA[[8]],listA[[9]],listA[[10]],listA[[11]])
listA[[10]]


# Relation between minimum distance and interaction type (all strains sum)
my_comparisons <- list(c("-4","-3"),c("-3","-2"),c("-2","-1"),
                       c("-1","0"),c("0","1"),c("1","2"),c("2","3"),c("3","4"),
                       c("-4","-2"),c("-3","-1"),c("-2","0"),c("-1","1"),
                       c("0","2"),c("1","3"),c("2","4"))
  
net <- df_tot$network %>% unique()
listA <- list()
for (i in 1:11){
  font_size <- 2
  listA[[i]] <- df_tot %>% filter(network==net[i]) %>% 
    ggboxplot(x = "sum_g", y = "path.length",
              fill = "sum_g")+ 
    stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = -3,size=font_size) +
    stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
    theme_minimal()+
    scale_fill_brewer(palette="Set2")+
    xlab("Type of interaction")+
    ylab("Minimum path length between proteins")+
    theme(legend.position = "none")+
    ggtitle(net[i])
  
  listA[[i]]$layers[[2]]$aes_params$textsize <- font_size
  
}
wrap_plots(listA[[1]],listA[[2]],listA[[3]],listA[[4]],listA[[5]],listA[[6]],
           listA[[7]],listA[[8]],listA[[9]],listA[[10]],listA[[11]])
listA[[10]]



# Relation between connectivity and interaction type (ebw). 
net <- df_tot$network %>% unique()
listA <- list()
for (i in 1:11){
  listA[[i]] <- df_tot %>% filter(network==net[i]) %>% 
    ggboxplot(x = "int_sign_ebw", y = "K.edge",
              fill = "int_sign_ebw")+ 
    stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = -3) +
    stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
    theme_minimal()+
    scale_fill_brewer(palette="Set2")+
    xlab("Type of interaction")+
    ylab("K edge connectivity between proteins")+
    theme(legend.position = "none")+
    ggtitle(net[i])
}

wrap_plots(listA[[1]],listA[[2]],listA[[3]],listA[[4]],listA[[5]],listA[[6]],
           listA[[7]],listA[[8]],listA[[9]],listA[[10]],listA[[11]])
listA[[10]]


# Relation between connectivity and interaction type (all strains sum)
my_comparisons <- list(c("-4","-3"),c("-3","-2"),c("-2","-1"),
                       c("-1","0"),c("0","1"),c("1","2"),c("2","3"),c("3","4"),
                       c("-4","-2"),c("-3","-1"),c("-2","0"),c("-1","1"),
                       c("0","2"),c("1","3"),c("2","4"))

net <- df_tot$network %>% unique()
listA <- list()
for (i in 1:11){
  font_size <- 2
  listA[[i]] <- df_tot %>% filter(network==net[i]) %>% 
    ggboxplot(x = "sum_g", y = "K.edge",
              fill = "sum_g")+ 
    stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = -3,size=font_size) +
    stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
    theme_minimal()+
    scale_fill_brewer(palette="Set2")+
    xlab("Type of interaction")+
    ylab("K edge connectivity between proteins")+
    theme(legend.position = "none")+
    ggtitle(net[i])
  
  listA[[i]]$layers[[2]]$aes_params$textsize <- font_size
  
}
wrap_plots(listA[[1]],listA[[2]],listA[[3]],listA[[4]],listA[[5]],listA[[6]],
           listA[[7]],listA[[8]],listA[[9]],listA[[10]],listA[[11]])
listA[[10]]


# Relation between node degree and interaction type (ebw). 
net <- df_tot$network %>% unique()
listA <- list()
for (i in 1:11){
  listA[[i]] <- df_tot %>% filter(network==net[i]) %>% 
    ggboxplot(x = "int_sign_ebw", y = "mean_deg",
              fill = "int_sign_ebw")+ 
    stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = -3) +
    stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
    theme_minimal()+
    scale_fill_brewer(palette="Set2")+
    xlab("Type of interaction")+
    ylab("Mean degree between proteins")+
    theme(legend.position = "none")+
    ggtitle(net[i])
}

wrap_plots(listA[[1]],listA[[2]],listA[[3]],listA[[4]],listA[[5]],listA[[6]],
           listA[[7]],listA[[8]],listA[[9]],listA[[10]],listA[[11]])
listA[[10]]

# Relation between node degree and interaction type (all strains sum)
my_comparisons <- list(c("-4","-3"),c("-3","-2"),c("-2","-1"),
                       c("-1","0"),c("0","1"),c("1","2"),c("2","3"),c("3","4"),
                       c("-4","-2"),c("-3","-1"),c("-2","0"),c("-1","1"),
                       c("0","2"),c("1","3"),c("2","4"))

net <- df_tot$network %>% unique()
listA <- list()
for (i in 1:11){
  font_size <- 2
  listA[[i]] <- df_tot %>% filter(network==net[i]) %>% 
    ggboxplot(x = "sum_g", y = "mean_deg",
              fill = "sum_g")+ 
    stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = -3,size=font_size) +
    stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
    theme_minimal()+
    scale_fill_brewer(palette="Set2")+
    xlab("Type of interaction")+
    ylab("Mean degree between proteins")+
    theme(legend.position = "none")+
    ggtitle(net[i])
  
  listA[[i]]$layers[[2]]$aes_params$textsize <- font_size
  
}
wrap_plots(listA[[1]],listA[[2]],listA[[3]],listA[[4]],listA[[5]],listA[[6]],
           listA[[7]],listA[[8]],listA[[9]],listA[[10]],listA[[11]])
listA[[10]]

# % Connected ~ Type
d1 <- df_tot %>%
  filter(network=="Small/medium-scale PPI") %>% 
  select(drug_pair,U_ID,connection_onoff,int_sign_ebw,sum_g)  %>% 
  mutate(connect_01=ifelse(connection_onoff=="disconnected",0,
                      ifelse(connection_onoff=="connected",1,NA))) %>% 
  mutate(connect_01=as.numeric(connect_01)) %>% 
  group_by(int_sign_ebw) %>% 
  summarise(mean_conn=mean(connect_01)) %>% as.data.frame()
  
d1 %>% ggplot(aes(x = int_sign_ebw, y = mean_conn))+
  geom_point()+
  theme_minimal()+
  scale_fill_brewer(palette="Set2")+
  xlab("Type of interaction")+
  ylab("Percentage connection")+
  theme(legend.position = "none")
  
# int_sign_ebw mean_conn
# 1      Synergy 0.6631579
# 2   Additivity 0.8236301
# 3   Antagonism 1.0000000

############################################
# For each cluster, calculate:
# number of unique DDIs per cluster
df_tot %>%
  filter(network=="Small/medium-scale PPI") %>% 
  select(sigma.rate,clusters,drug_pair,connection_onoff,int_sign_ebw,sum_g) %>% 
  distinct() %>% 
  group_by(sigma.rate) %>%
  summarise(count = n_distinct(drug_pair)) %>% 
  ggplot(aes(x=sigma.rate,y=count,group=sigma.rate))+geom_point()+theme_minimal()
# number of unique target pairs
df_tot %>%
  filter(network=="Small/medium-scale PPI") %>% 
  select(sigma.rate,clusters,drug_pair,U_ID,connection_onoff,int_sign_ebw,sum_g) %>% 
  group_by(sigma.rate) %>%
  summarise(count = n_distinct(U_ID)) %>% 
  ggplot(aes(x=sigma.rate,y=count))+geom_point()+theme_minimal()

# average min path length
df_tot %>%
  filter(network=="Small/medium-scale PPI") %>% 
  select(sigma.rate,clusters,drug_pair,connection_onoff,int_sign_ebw,sum_g,
         path.length)  %>% 
  distinct() %>% 
  ggplot(aes(x=sigma.rate,y=path.length,group=sigma.rate))+geom_boxplot()+
  theme_minimal()

# average k edge
df_tot %>%
  filter(network=="Small/medium-scale PPI") %>% 
  select(sigma.rate,clusters,drug_pair,connection_onoff,int_sign_ebw,sum_g,
         K.edge)  %>% 
  distinct() %>% 
  ggplot(aes(x=sigma.rate,y=K.edge,group=sigma.rate))+geom_boxplot()+
  geom_jitter()+theme_minimal()
# average node degree
df_tot %>%
  filter(network=="Small/medium-scale PPI") %>% 
  select(sigma.rate,clusters,drug_pair,connection_onoff,int_sign_ebw,sum_g,
         mean_deg)  %>% 
  distinct() %>% 
  ggplot(aes(x=sigma.rate,y=mean_deg,group=sigma.rate))+geom_boxplot()+
  geom_jitter()+theme_minimal()
# % of connected vs disconnected

# % synergies, % additivities, % antagonisms

# average sum_g
df_tot %>%
  filter(network=="Small/medium-scale PPI") %>% 
  select(sigma.rate,clusters,drug_pair,connection_onoff,int_sign_ebw,sum_g,
         sum_g)  %>% 
  distinct() %>% 
  ggplot(aes(x=sigma.rate,y=sum_g,group=sigma.rate))+geom_boxplot()+
  geom_jitter()+theme_minimal()


# mean path by type of interaction (Syn,Ant,Add,Add-Ant,Ant-Syn,Add-Syn,Add-Ant-Syn)
r <- c()
for (i in 1:dim(df_tot)[1]){
  x <- df_tot[i,35:40] %>% t() %>% as.vector()
  x <- x %>% unique() %>% sort()
  r[i] <- paste0(x,collapse = "-")
}
df_tot$type <- r


# Mean.path
df_tot %>%
  filter(network=="Small/medium-scale PPI") %>% 
  select(sigma.rate,clusters,drug_pair,connection_onoff,type,path.length)  %>% 
  distinct() %>% 
  group_by(type) %>% 
  summarise(Mean.path=mean(path.length),
            sem.path=std.error(path.length)) %>% 
  ungroup() %>% 
  ggplot(aes(x=type,y=Mean.path,fill=type))+
  geom_bar(position=position_dodge(), stat="identity",
           colour='black',width=0.5)+
  geom_errorbar(position=position_dodge(0.5),aes(ymin=Mean.path-sem.path, ymax=Mean.path+sem.path), width=.2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Interaction type across all strains")+
  ylab("Minimum path distance")+
  theme(legend.position = "none")

# Mean.rate
df_tot %>%
  filter(network=="Small/medium-scale PPI") %>% 
  select(sigma.rate,clusters,drug_pair,connection_onoff,type,path.length)  %>% 
  distinct() %>% 
  group_by(type) %>% 
  summarise(Mean.rate=mean(sigma.rate),
            sem.rate=std.error(sigma.rate)) %>% 
  ungroup() %>% 
  ggplot(aes(x=type,y=Mean.rate,fill=type))+
  geom_bar(position=position_dodge(), stat="identity",
           colour='black',width=0.5)+
  geom_errorbar(position=position_dodge(0.5),aes(ymin=Mean.rate-sem.rate, ymax=Mean.rate+sem.rate), width=.2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Interaction type across all strains")+
  ylab("Sigma rate")+
  theme(legend.position = "none")


df_tot %>%
  filter(network=="Small/medium-scale PPI") %>% 
  select(sigma.rate,clusters,drug_pair,connection_onoff,type,path.length,K.edge)  %>% 
  distinct() %>% 
  group_by(type) %>% 
  summarise(Mean.K.edge=mean(K.edge,na.rm=T),
            sem.K.edge=std.error(K.edge)) %>% 
  ungroup() %>% 
  ggplot(aes(x=type,y=Mean.K.edge,fill=type))+
  geom_bar(position=position_dodge(), stat="identity",
           colour='black',width=0.5)+
  geom_errorbar(position=position_dodge(0.5),aes(ymin=Mean.K.edge-sem.K.edge, ymax=Mean.K.edge+sem.K.edge), width=.2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


df_tot %>%
  filter(network=="Small/medium-scale PPI") %>% 
  select(sigma.rate,clusters,drug_pair,connection_onoff,type,path.length,
         mean_deg)  %>% 
  distinct() %>% 
  group_by(type) %>% 
  summarise(Mean.mean_deg=mean(mean_deg,na.rm=T),
            sem.mean_deg=std.error(mean_deg)) %>% 
  ungroup() %>% 
  ggplot(aes(x=type,y=Mean.mean_deg,fill=type))+
  geom_bar(position=position_dodge(), stat="identity",
           colour='black',width=0.5)+
  geom_errorbar(position=position_dodge(0.5),aes(ymin=Mean.mean_deg-sem.mean_deg, ymax=Mean.mean_deg+sem.mean_deg), width=.2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



df_tot %>%
  filter(network=="Small/medium-scale PPI") %>% 
  select(sigma.rate,clusters,drug_pair,connection_onoff,type,path.length)  %>% 
  distinct() %>% 
  ggplot(aes(x=type,y=sigma.rate))+geom_boxplot()+
  geom_jitter(aes(col=clusters))+theme_minimal()

df_tot %>%
  filter(network=="Small/medium-scale PPI") %>% 
  select(sigma.rate,clusters,drug_pair,connection_onoff,type,path.length)  %>% 
  distinct() %>% 
  ggplot(aes(x=type,y=as.numeric(as.character(clusters))))+geom_boxplot()+
  geom_jitter(aes(col=clusters))+theme_minimal()


# Evolutionary rates ~ int_type/distance/connectivity/node [Circular]

# Relation between sigma rate and interaction type (ebw). 
net <- df_tot$network %>% unique()
listA <- list()
for (i in 1:11){
  font_size <- 3
  
  listA[[i]] <- df_tot %>% filter(network==net[i]) %>% 
    ggboxplot(x = "int_sign_ebw", y = "sigma.rate",
              fill = "int_sign_ebw")+ 
    stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = -3,size=font_size) +
    stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
    theme_minimal()+
    scale_fill_brewer(palette="Set2")+
    xlab("Type of interaction")+
    ylab("Sigma rate")+
    theme(legend.position = "none")+
    ggtitle(net[i])
  
  listA[[i]]$layers[[2]]$aes_params$textsize <- font_size
}
wrap_plots(listA[[1]],listA[[2]],listA[[3]],listA[[4]],listA[[5]],listA[[6]],
           listA[[7]],listA[[8]],listA[[9]],listA[[10]],listA[[11]])
listA[[10]]


# Relation between sigma rate and interaction type (all strains sum)
my_comparisons <- list(c("-4","-3"),c("-3","-2"),c("-2","-1"),
                       c("-1","0"),c("0","1"),c("1","2"),c("2","3"),c("3","4"),
                       c("-4","-2"),c("-3","-1"),c("-2","0"),c("-1","1"),
                       c("0","2"),c("1","3"),c("2","4"))

net <- df_tot$network %>% unique()
listA <- list()
for (i in 1:11){
  font_size <- 2
  listA[[i]] <- df_tot %>% filter(network==net[i]) %>% 
    ggboxplot(x = "sum_g", y = "sigma.rate",
              fill = "sum_g")+ 
    stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = -3,size=font_size) +
    stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
    theme_minimal()+
    scale_fill_brewer(palette="Set2")+
    xlab("Type of interaction")+
    ylab("Sigma rate")+
    theme(legend.position = "none")+
    ggtitle(net[i])
  
  listA[[i]]$layers[[2]]$aes_params$textsize <- font_size
  
}
wrap_plots(listA[[1]],listA[[2]],listA[[3]],listA[[4]],listA[[5]],listA[[6]],
           listA[[7]],listA[[8]],listA[[9]],listA[[10]],listA[[11]])
listA[[10]]

# Rates ~ Interaction Type (ebw) (all networks).
df_tot %>% ggline(x = "int_sign_ebw", y = "sigma.rate", add = "mean_se",
       group = "network", palette = "jco")+
  stat_compare_means(aes(group = network), label = "p.signif", 
                     label.y = c(40, 40, 40))+
  facet_wrap(~network)+
  xlab("Type of interaction")+
  ylab("Sigma rate")+
  theme_minimal()

# Rates ~ Interaction Type (all strains) (all networks).
df_tot %>% ggline(x = "sum_g", y = "sigma.rate", add = "mean_se",
                  group = "network", palette = "jco")+
  stat_compare_means(aes(group = network), label = "p.signif", 
                     label.y = c(40, 40, 40))+
  facet_wrap(~network)+
  xlab("Type of interaction")+
  ylab("Sigma rate")+
  theme_minimal()

# Rates ~ Interaction Type (ebw) (PPI). [This is circular]
df_tot %>%
  filter(network=="Small/medium-scale PPI") %>% 
  ggline(x = "int_sign_ebw", y = "sigma.rate", add = "mean_se",
                  group = "network", palette = "jco")+
  stat_compare_means(aes(group = network), label = "p.signif", 
                     label.y = c(40, 40, 40))+
  facet_wrap(~network)+
  xlab("Type of interaction")+
  ylab("Sigma rate")+
  theme_minimal()

# Rates ~ Interaction Type (all strains) (PPI). [This is circular]
df_tot %>%
  filter(network=="Small/medium-scale PPI") %>% 
  ggline(x = "sum_g", y = "sigma.rate", add = "mean_se",
         group = "network", palette = "jco")+
  stat_compare_means(aes(group = network), label = "p.signif", 
                     label.y = c(40, 40, 40))+
  facet_wrap(~network)+
  xlab("Type of interaction")+
  ylab("Sigma rate")+
  theme_minimal()

# Rates ~ Distance (PPI).
df_tot %>% 
  filter(network=="Small/medium-scale PPI") %>% 
  mutate(path.length.corr=ifelse(path.length==0,"0",
                                 ifelse(path.length==1,"1",
                                        ifelse(path.length==2,"2",
                                               ifelse(path.length==3,"3",
                                                      ifelse(path.length==4,"4",
                                                             ifelse(path.length==5,"5","6+"))))))) %>% 
  mutate(path.length.corr=factor(path.length.corr,levels=c("0","1","2","3","4","5","6+"))) %>% 
  ggline(x = "path.length.corr", y = "sigma.rate", add = "mean_se", 
         group = "network",palette = "jco")+
  stat_compare_means(aes(group = network), label = "p.signif", 
                     label.y = c(40, 40, 40))+
  facet_wrap(~network)+
  xlab("Minimum distance between targets")+
  ylab("Sigma rate")+
  theme_minimal()


# Rates ~ Distance (all networks).
df_tot %>%  mutate(path.length.corr=ifelse(path.length==0,"0",
                                 ifelse(path.length==1,"1",
                                        ifelse(path.length==2,"2",
                                               ifelse(path.length==3,"3",
                                                      ifelse(path.length==4,"4",
                                                             ifelse(path.length==5,"5","6+"))))))) %>% 
  mutate(path.length.corr=factor(path.length.corr,levels=c("0","1","2","3","4","5","6+"))) %>% 
  ggline(x = "path.length.corr", y = "sigma.rate", add = "mean_se",
                  group = "network", palette = "jco")+
  stat_compare_means(aes(group = network), label = "p.signif", 
                     label.y = c(40, 40, 40))+
  facet_wrap(~network)+
  xlab("Minimum distance between targets")+
  ylab("Sigma rate")+
  theme_minimal()

# Rates ~ Connectivity (all networks).
df_tot %>% ggline(x = "K.edge", y = "sigma.rate", add = "mean_se",
                  group = "network", palette = "jco")+
  stat_compare_means(aes(group = network), label = "p.signif", 
                     label.y = c(40, 40, 40))+
  facet_wrap(~network,scales="free")+
  xlab("K-edge connectivity between proteins")+
  ylab("Sigma rate")+
  theme_minimal()

# Rates ~ Node degree (all networks).
df_tot %>% ggline(x = "mean_deg", y = "sigma.rate", add = "mean_se",
                  group = "network", palette = "jco")+
  stat_compare_means(aes(group = network), label = "p.signif", 
                     label.y = c(40, 40, 40))+
  facet_wrap(~network)+
  xlab("Average node degree between target proteins")+
  ylab("Sigma rate")+
  theme_minimal()

# Rate ~ Connectivity (all networks)
df_tot %>% ggline(x = "connection_onoff", y = "sigma.rate", add = "mean_se",
                  col = "network")+
  stat_compare_means(aes(group = network), label = "p.signif", 
                     label.y = c(40, 40, 40))+
  xlab("Interconnectedness of nodes")+
  ylab("Sigma rate")+
  theme_minimal()

# Rate ~ Adjacency (all networks)
df_tot %>% ggline(x = "adjancency", y = "sigma.rate", add = "mean_se",
                  col = "network")+
  stat_compare_means(aes(group = network), label = "p.signif", 
                     label.y = c(40, 40, 40))+
  xlab("Adjacency of nodes")+
  ylab("Sigma rate")+
  theme_minimal()

# Rate ~ interaction type (across species)
# sum_g just adds the scores of interaction types (-1:antagonism,0:additive,+1:synergy)
df_tot %>%
  filter(network=="Small/medium-scale PPI") %>% 
  ggline(x = "sum_g", y = "sigma.rate", add = "mean_se",
         group = "network", palette = "jco")+
  stat_compare_means(aes(group = network), label = "p.signif", 
                     label.y = c(40, 40, 40))+
  facet_wrap(~network)+
  xlab("Type of interaction")+
  ylab("Sigma rate")+
  theme_minimal()

df_tot %>%
  filter(network=="Small/medium-scale PPI") %>% 
  select(sum_g) %>% pull() %>% hist()

# Cluster ~ interaction type
df_tot %>%
  filter(network=="Small/medium-scale PPI") %>% 
  ggplot(aes(x=clusters,y=sd_g,#sd_g, sum_g
             group=clusters,
             fill=sigma.rate))+
  geom_boxplot()+
  theme_minimal()


# clusters with low sigma rates have different categories.
# clusters with high sigma rates come from the same category.
# Cluster 13 has most synergies, cluster 2 has most antagonisms. They have the highest
# rates.
df_tot %>% 
  select(int_sign_ebw,int_sign_ecr,
         int_sign_seo,int_sign_stm,
         int_sign_pae,int_sign_pau,
         sigma.rate) %>% 
  pivot_longer(names_to="species",values_to = "signs",1:6) %>% 
  ggplot(aes(x=signs,y=sigma.rate))+
  geom_point()+
  geom_violin()

