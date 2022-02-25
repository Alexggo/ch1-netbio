#' ---
#' title: "5.2.Do connected proteins have faster or slower rates?"
#' author: "Alex Gil"
#' 
library(pacman)
p_load(igraph,tidyverse,BioNet,ape,geomorph,matrixStats,
       ComplexHeatmap,ggpubr,plotrix,tidymodels)


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
              "LC.INT.EcoliNet.764gene.1073link.txt","PG.INT.EcoliNet.v1.1817gene.17504link.txt",
              "MODEL1108160000_edgelist.txt"
)



rates <- "dataset1_ratios.csv"
rates <- read.csv(file.path("data/4.PhylogeneticComparativeMethods",rates))

list1 <- list()
for (i in 1:12){
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
    pivot_longer(names_to="Uniprot2",values_to="adjacency",2:(dim(distnew2)[2]))%>%
    filter(!is.na(adjacency)) %>% 
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
names(list1) <- filename


# List 1 contains the complete dataset for each network.
# df_tot contains all networks combined.
df_tot <- bind_rows(list1, .id = "network") %>% 
  mutate(adjacency=as.factor(adjacency)) %>% 
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
                                                                                        ifelse(network=="PG.INT.EcoliNet.v1.1817gene.17504link.txt","Similar phylogenetic profiles",
                                                                                               ifelse(network=="MODEL1108160000_edgelist.txt","Metabolic model","NA"))))))))))))) %>% 
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
net <- df_tot$network %>% unique()

# Minimum path ~ interaction type (ebw). 
my_comparisons <- list( c("Synergy", "Additivity"), c("Additivity", "Antagonism"), c("Synergy", "Antagonism") )
net <- df_tot$network %>% unique()
listA <- list()
for (i in 1:12){
  
  font_size <- 2
  listA[[i]] <- df_tot %>% filter(network==net[i]) %>% 
    ggplot(aes(x=int_sign_ebw,y=path.length,fill=int_sign_ebw))+
    geom_violin()+
    geom_boxplot(width=0.01,fill="white")+
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
           listA[[7]],listA[[8]],listA[[9]],listA[[10]],listA[[11]],listA[[12]])
listA[[10]]
listA[[12]]
wrap_plots(listA[[10]],listA[[12]])
# Synergies occur closer in the network.


# Connectivity ~ interaction type (ebw).
my_comparisons <- list( c("Synergy", "Additivity"), c("Additivity", "Antagonism"), c("Synergy", "Antagonism") )
net <- df_tot$network %>% unique()
listA <- list()
for (i in 1:12){
  listA[[i]] <- df_tot %>% filter(network==net[i]) %>% 
    ggplot(aes(x=int_sign_ebw,y=K.edge,fill=int_sign_ebw))+
    geom_violin()+
    geom_boxplot(width=0.01,fill="white")+
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
           listA[[7]],listA[[8]],listA[[9]],listA[[10]],listA[[11]],listA[[12]])
listA[[10]]
listA[[12]]
wrap_plots(listA[[10]],listA[[12]])
#Con is lower in antagonisms but only in metabolic model.


# Node degree ~ interaction type (ebw). 
my_comparisons <- list( c("Synergy", "Additivity"), c("Additivity", "Antagonism"), c("Synergy", "Antagonism") )
net <- df_tot$network %>% unique()
listA <- list()
for (i in 1:12){
  listA[[i]] <- df_tot %>% filter(network==net[i]) %>% 
    ggplot(aes(x = int_sign_ebw, y = mean_deg,
              fill = int_sign_ebw))+ 
    geom_violin()+
    geom_boxplot(width=0.01,fill="white")+
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
           listA[[7]],listA[[8]],listA[[9]],listA[[10]],listA[[11]],listA[[12]])
listA[[10]]
listA[[12]]
wrap_plots(listA[[10]],listA[[12]])
# Significant differences (mean degree) between interactions in the metabolic network.



# % Connected ~ Type
# For each type of DDI, how many have its targets connected?
listA <- list()
listB <- list()
for (i in 1:12){
d1 <- df_tot %>%
    filter(network==net[i]) %>% 
    select(drug_pair,U_ID,connection_onoff,int_sign_ebw,sum_g)  %>% 
    mutate(connect_01=ifelse(connection_onoff=="disconnected",0,
                             ifelse(connection_onoff=="connected",1,NA))) %>% 
    mutate(connect_01=as.numeric(connect_01)) %>% 
    group_by(int_sign_ebw) 
  listA[[i]] <- d1
  
d2 <- d1 %>% ggplot(aes(x = int_sign_ebw, fill = connection_onoff))+
    geom_bar(position = "dodge")+
    theme_minimal()+
    scale_fill_brewer(palette="Set2")+
    xlab("Type of interaction")+
    ylab("Number of target connections")+
    theme(legend.position = "none")+
    ggtitle(net[i])
  listB[[i]] <- d2
}
wrap_plots(listB[[1]],listB[[2]],listB[[3]],listB[[4]],listB[[5]],listB[[6]],
           listB[[7]],listB[[8]],listB[[9]],listB[[10]],listB[[11]],listB[[12]])

wrap_plots(listB[[10]],listB[[12]])
#connected green, disconnected orange.
#Connected is more common than disconnected in all categories.

# % Adjacency ~ Type
# For each type of DDI, how many have its targets connected?
listA <- list()
listB <- list()
for (i in 1:12){
  d1 <- df_tot %>%
    filter(network==net[i]) %>% 
    select(drug_pair,U_ID,adjacency,int_sign_ebw,sum_g)  %>% 
    mutate(adja=ifelse(adjacency==0,"not adjacent",
                             ifelse(adjacency==1,"adjacent",NA))) %>% 
    group_by(int_sign_ebw) 
  listA[[i]] <- d1
  
  d2 <- d1 %>% ggplot(aes(x = int_sign_ebw, fill = adja))+
    geom_bar(position = "dodge")+
    theme_minimal()+
    scale_fill_brewer(palette="Set2")+
    xlab("Type of interaction")+
    ylab("Number of target connections")+
    theme(legend.position = "bottom")+
    ggtitle(net[i])
  listB[[i]] <- d2
}
wrap_plots(listB[[1]],listB[[2]],listB[[3]],listB[[4]],listB[[5]],listB[[6]],
           listB[[7]],listB[[8]],listB[[9]],listB[[10]],listB[[11]],listB[[12]])

wrap_plots(listB[[10]],listB[[12]])
#Adjacent green, not adjacent orange
#Not adjacent is more common than adjacent, in all categories.

# Rate ~ Conn/Adj
my_comparisons <- list(c("0","1"))
listCon <- list()
listAdj <- list()
for (i in 1:12){
  # Rate ~ connected/non connected
  listCon[[i]] <- df_tot %>%
  filter(network==net[i]) %>% 
  ggplot(aes(y=sigma.rate,x=adjacency,
             fill=adjacency))+
  geom_boxplot()+
    stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = -3,size=font_size) +
    stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
  theme_minimal()+
  theme(legend.position = "none")+
  ggtitle(net[i])+
  ylab("Sigma rate")+
  xlab("Connectedness")
# Rates ~ adjacent/non adjacent
  listAdj[[i]] <- df_tot %>%
  filter(network==net[i]) %>% 
  ggplot(aes(y=sigma.rate,x=connection_onoff,
             fill=connection_onoff))+
  geom_boxplot()+
    stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = -3,size=font_size) +
    stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
  theme_minimal()+
  theme(legend.position = "none")+
  ggtitle(net[i])+
  ylab("Sigma rate")+
  xlab("Adjacency")
}

#Sup. Fig 3,4
wrap_plots(listCon[[10]],listCon[[12]])
wrap_plots(listAdj[[10]],listAdj[[12]])





############################################
# For each cluster, calculate:

# number of unique DDIs per cluster
df_tot  %>% 
  select(network,sigma.rate,clusters,drug_pair,connection_onoff,int_sign_ebw,sum_g) %>% 
  distinct() %>% 
  group_by(sigma.rate,network) %>%
  summarise(count = n_distinct(drug_pair)) %>% 
  ggplot(aes(x=sigma.rate,y=count,group=sigma.rate))+geom_point()+theme_minimal()+
  facet_wrap(.~network)+
  ylab("Number of DDIs")
# number of unique target pairs
df_tot %>% 
  select(network,sigma.rate,clusters,drug_pair,U_ID,connection_onoff,int_sign_ebw,sum_g) %>% 
  group_by(network,sigma.rate) %>%
  summarise(count = n_distinct(U_ID)) %>% 
  ggplot(aes(x=sigma.rate,y=count))+geom_point()+
  facet_wrap(~network)+
  theme_minimal()+
  ylab("Number of targets")

# average min path length per cluster
df_tot %>%
  filter(network=="Small/medium-scale PPI") %>% 
  select(sigma.rate,clusters,drug_pair,connection_onoff,int_sign_ebw,sum_g,
         path.length)  %>% 
  distinct() %>% 
  ggplot(aes(x=sigma.rate,y=path.length,group=sigma.rate))+geom_boxplot()+
  theme_minimal()

# average k edge per cluster
df_tot %>%
  filter(network=="Small/medium-scale PPI") %>% 
  select(sigma.rate,clusters,drug_pair,connection_onoff,int_sign_ebw,sum_g,
         K.edge)  %>% 
  distinct() %>% 
  ggplot(aes(x=sigma.rate,y=K.edge,group=sigma.rate))+
  geom_boxplot()+
  theme_minimal()
# average node degree
df_tot %>%
  filter(network=="Small/medium-scale PPI") %>% 
  select(sigma.rate,clusters,drug_pair,connection_onoff,int_sign_ebw,sum_g,
         mean_deg)  %>% 
  distinct() %>% 
  ggplot(aes(x=sigma.rate,y=mean_deg,group=sigma.rate))+geom_boxplot()+
  theme_minimal()

# average sum_g
df_tot %>%
  filter(network=="Small/medium-scale PPI") %>% 
  select(sigma.rate,clusters,drug_pair,connection_onoff,int_sign_ebw,sum_g,
         sum_g)  %>% 
  distinct() %>% 
  ggplot(aes(x=sigma.rate,y=sum_g,group=sigma.rate))+geom_boxplot()+
  geom_jitter()+theme_minimal()
# clusters with higher rates are mainly a single int type.


# mean path ~ interaction (Syn,Ant,Add,Add-Ant,Ant-Syn,Add-Syn,Add-Ant-Syn)
r <- c()
for (i in 1:dim(df_tot)[1]){
  x <- df_tot[i,33:38] %>% t() %>% as.vector()
  x <- x %>% unique() %>% sort()
  r[i] <- paste0(x,collapse = "-")
}
df_tot$type <- r


# Mean.path
df_tot %>%
  filter(network==net[10]|network==net[12]) %>% 
  select(network,sigma.rate,clusters,drug_pair,connection_onoff,type,path.length)  %>% 
  distinct() %>% 
  group_by(type,network) %>% 
  summarise(Mean.path=mean(path.length),
            sem.path=std.error(path.length)) %>% 
  ungroup() %>% 
  ggplot(aes(x=type,y=Mean.path,fill=network,group=network))+
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
  filter(network==net[10]|network==net[12]) %>% 
  select(network,sigma.rate,clusters,drug_pair,connection_onoff,type,path.length)  %>% 
  distinct() %>% 
  group_by(network,type) %>% 
  summarise(Mean.rate=mean(sigma.rate),
            sem.rate=std.error(sigma.rate)) %>% 
  ungroup() %>% 
  ggplot(aes(x=type,y=Mean.rate,fill=network))+
  geom_bar(position=position_dodge(), stat="identity",
           colour='black',width=0.5)+
  geom_errorbar(position=position_dodge(0.5),aes(ymin=Mean.rate-sem.rate, ymax=Mean.rate+sem.rate), width=.2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Interaction type across all strains")+
  ylab("Sigma rate")+
  theme(legend.position = "left")


df_tot %>%
  filter(network==net[10]|network==net[12]) %>% 
  select(network,sigma.rate,clusters,drug_pair,connection_onoff,type,path.length,K.edge)  %>% 
  distinct() %>% 
  group_by(type,network) %>% 
  summarise(Mean.K.edge=mean(K.edge,na.rm=T),
            sem.K.edge=std.error(K.edge)) %>% 
  ungroup() %>% 
  ggplot(aes(x=type,y=Mean.K.edge,fill=network))+
  geom_bar(position=position_dodge(), stat="identity",
           colour='black',width=0.5)+
  geom_errorbar(position=position_dodge(0.5),aes(ymin=Mean.K.edge-sem.K.edge, ymax=Mean.K.edge+sem.K.edge), width=.2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


df_tot %>%
  filter(network==net[10]|network==net[12]) %>% 
  select(network,sigma.rate,clusters,drug_pair,connection_onoff,type,path.length,
         mean_deg)  %>% 
  distinct() %>% 
  group_by(type,network) %>% 
  summarise(Mean.mean_deg=mean(mean_deg,na.rm=T),
            sem.mean_deg=std.error(mean_deg)) %>% 
  ungroup() %>% 
  ggplot(aes(x=type,y=Mean.mean_deg,fill=network))+
  geom_bar(position=position_dodge(), stat="identity",
           colour='black',width=0.5)+
  geom_errorbar(position=position_dodge(0.5),aes(ymin=Mean.mean_deg-sem.mean_deg, ymax=Mean.mean_deg+sem.mean_deg), width=.2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Evolutionary rates ~ int_type/distance/connectivity/node 

# Relation between sigma rate and interaction type (ebw). 
net <- df_tot$network %>% unique()
listA <- list()
for (i in 1:12){
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
           listA[[7]],listA[[8]],listA[[9]],listA[[10]],listA[[11]],listA[[12]])
listA[[10]]
wrap_plots(listA[[10]],listA[[12]])

# Relation between sigma rate and interaction type (all strains sum)
my_comparisons <- list(c("-4","-3"),c("-3","-2"),c("-2","-1"),
                       c("-1","0"),c("0","1"),c("1","2"),c("2","3"),c("3","4"),
                       c("-4","-2"),c("-3","-1"),c("-2","0"),c("-1","1"),
                       c("0","2"),c("1","3"),c("2","4"))

net <- df_tot$network %>% unique()
listA <- list()
for (i in 1:12){
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
           listA[[7]],listA[[8]],listA[[9]],listA[[10]],listA[[11]],listA[[12]])
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

# Rates ~ Interaction Type (ebw) (PPI).
df_tot %>%
  filter(network==net[10]|network==net[12]) %>% 
  ggline(x = "int_sign_ebw", y = "sigma.rate", add = "mean_se",
                  group = "network", col="network")+
  xlab("Type of interaction")+
  ylab("Sigma rate")+
  theme_minimal()+
  theme(legend.position = "none")

# Rates ~ Interaction Type (all strains) (PPI).
df_tot %>%
  filter(network==net[10]|network==net[12]) %>% 
  ggline(x = "sum_g", y = "sigma.rate", add = "mean_se",
         group = "network",col="network")+
  xlab("Type of interaction")+
  ylab("Sigma rate")+
  theme_minimal()+
  theme(legend.position = "none")

# Rates ~ Distance (PPI).
df_tot %>% 
  filter(network==net[10]|network==net[12]) %>% 
  mutate(path.length.corr=ifelse(path.length==0,"0",
                                 ifelse(path.length==1,"1",
                                        ifelse(path.length==2,"2",
                                               ifelse(path.length==3,"3",
                                                      ifelse(path.length==4,"4",
                                                             ifelse(path.length==5,"5","6+"))))))) %>% 
  mutate(path.length.corr=factor(path.length.corr,levels=c("0","1","2","3","4","5","6+"))) %>% 
  ggline(x = "path.length.corr", y = "sigma.rate", add = "mean_se", 
         group = "network",col="network")+
  xlab("Minimum distance between targets")+
  ylab("Sigma rate")+
  theme_minimal()+
  theme(legend.position = "none")


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
df_tot %>% ggline(x = "adjacency", y = "sigma.rate", add = "mean_se",
                  col = "network")+
  stat_compare_means(aes(group = network), label = "p.signif", 
                     label.y = c(40, 40, 40))+
  xlab("Adjacency of nodes")+
  ylab("Sigma rate")+
  theme_minimal()

# Rate ~ interaction type (across species)
# sum_g just adds the scores of interaction types (-1:antagonism,0:additive,+1:synergy)
df_tot %>%
  filter(network==net[10]) %>% 
  ggline(x = "sum_g", y = "sigma.rate", add = "mean_se",
         group = "network", palette = "jco")+
  stat_compare_means(aes(group = network), label = "p.signif", 
                     label.y = c(40, 40, 40))+
  facet_wrap(~network)+
  xlab("Type of interaction")+
  ylab("Sigma rate")+
  theme_minimal()

