#' ---
#' title: "5.2.Do connected proteins have faster or slower rates?"
#' author: "Alex Gil"
#' 
library(igraph)
library(tidyverse)
library(BioNet)
library(ape)
library(geomorph)
library(broom)
library(patchwork)
library(ggridges)
library(ggpubr)

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
  t <- file.path("data/5.Targets_NetworkDistance",paste0(filename[i],"net.txt"))
  x1 <- read_csv(t)
  x1 <- x1 %>% select(Uni_N1,Uni_N2)
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

df_tot <- bind_rows(list1, .id = "id") %>% 
  mutate(adjancency=as.factor(adjancency)) %>% 
  mutate(connection_groups=cut_interval(K.edge, 6)) %>% 
  mutate(connection_groups=ifelse(is.na(connection_groups),0,connection_groups)) %>% 
  mutate(connection_onoff=ifelse(is.na(K.edge),"disconnected","connected")) %>% 
  rowwise() %>% 
  mutate(mean.ddi=mean(c(ebw,ecr,seo,stm,pae,pau))) %>% 
  mutate(sd.ddi=sd(c(ebw,ecr,seo,stm,pae,pau))) %>% 
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
  mutate(int_sign_ebw = factor(int_sign_ebw,levels = c("Synergy","Additivity","Antagonism"))) 
  

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

# Rates ~ Interaction Type.
df_tot %>% ggline(x = "int_sign_ebw", y = "sigma.rate", add = "mean_se",
       group = "network", palette = "jco")+
  stat_compare_means(aes(group = network), label = "p.signif", 
                     label.y = c(40, 40, 40))+
  facet_wrap(~network)+
  xlab("Type of interaction")+
  ylab("Sigma rate")+
  theme_minimal()

df_tot %>% ggline(x = "path.length", y = "sigma.rate", add = "mean_se",
                  group = "network", palette = "jco")+
  stat_compare_means(aes(group = network), label = "p.signif", 
                     label.y = c(40, 40, 40))+
  facet_wrap(~network)+
  xlab("Minimum distance between targets")+
  ylab("Sigma rate")+
  theme_minimal()

df_tot %>% ggline(x = "K.edge", y = "sigma.rate", add = "mean_se",
                  group = "network", palette = "jco")+
  stat_compare_means(aes(group = network), label = "p.signif", 
                     label.y = c(40, 40, 40))+
  facet_wrap(~network,scales="free")+
  xlab("K-edge connectivity between proteins")+
  ylab("Sigma rate")+
  theme_minimal()


df_tot %>% ggline(x = "mean_deg", y = "sigma.rate", add = "mean_se",
                  group = "network", palette = "jco")+
  stat_compare_means(aes(group = network), label = "p.signif", 
                     label.y = c(40, 40, 40))+
  facet_wrap(~network)+
  xlab("Average node degree between target proteins")+
  ylab("Sigma rate")+
  theme_minimal()

#Rate~Connectivity
df_tot %>% ggline(x = "connection_onoff", y = "sigma.rate", add = "mean_se",
                  group = "network", palette = "jco")+
  stat_compare_means(aes(group = network), label = "p.signif", 
                     label.y = c(40, 40, 40))+
  facet_wrap(~network)+
  xlab("Interconnectedness of nodes")+
  ylab("Sigma rate")+
  theme_minimal()

df_tot %>% ggplot(aes(x=clusters,fill=connection_onoff))+
  geom_bar(position = "fill")+
  theme_minimal()+
  facet_wrap(.~network)

df_tot %>% ggplot(aes(x=clusters,fill=connection_onoff))+
  geom_bar(position = "identity")+
  theme_minimal()+
  facet_wrap(.~network)

df_tot %>% ggplot(aes(x=connection_onoff,y=sigma.rate,
                      fill=connection_onoff))+
  geom_boxplot()+
  stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red") +
  theme_bw()+
  facet_wrap(.~network)


# Rates~Adjacency 
df_tot %>% ggline(x = "adjancency", y = "sigma.rate", add = "mean_se",
                  group = "network", palette = "jco")+
  stat_compare_means(aes(group = network), label = "p.signif", 
                     label.y = c(40, 40, 40))+
  facet_wrap(~network)+
  xlab("Adjacency")+
  ylab("Sigma rate")+
  theme_minimal()

df_tot %>% ggplot(aes(x=clusters,fill=adjancency))+
  geom_bar(position = "fill")+
  theme_minimal()+
  facet_wrap(.~network)

df_tot %>% ggplot(aes(x=clusters,fill=adjancency))+
  geom_bar(position = "identity")+
  theme_minimal()+
  facet_wrap(.~network)

df_tot %>% filter(int_sign_ebw=="Antagonism") %>% ggplot(aes(x=clusters,fill=int_sign_ebw))+
  geom_bar(position = "identity")+
  theme_minimal()+
  facet_wrap(.~network)







