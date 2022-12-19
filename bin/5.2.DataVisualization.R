library(pacman)
p_load(igraph,tidyverse,BioNet,ape,geomorph,matrixStats,
       ComplexHeatmap,ggpubr,plotrix,tidymodels,patchwork,here)
library(EnvStats) # for adding sample size

filenames <- c("EcoCyc.goldstandardset.csv", "EcoliNet.v1.csv",
               "GN.INT.EcoliNet.3568gene.23439link.csv",
               "GO-BP.goldstandardset.csv",
               "CC.EcoliNet.v1.2296gene.50528link.csv",
               "CX.INT.EcoliNet.v1.4039gene.67494link.csv",
               "DC.EcoliNet.2283gene.9643link.csv",
               "EcoliNet.v1.benchmark.csv",
               "HT.INT.EcoliNet.3209gene.15543link.csv",
               "LC.INT.EcoliNet.764gene.1073link.csv",
               "PG.INT.EcoliNet.v1.1817gene.17504link.csv")

net1 <- c("Co-functional(EcoCyc)", "EN:allnetworks",
         "Similar genomic context", "Co-functional (GO-BP)",
         "Co-citation", "Co-expression",
         "Co-occurence of prot. domains", "Co-functional (EcoCyc/GO-BP)",
         "High-throughput PPI", "Small/medium-scale PPI",
         "Similar phylogenetic profiles")

# VISUALIZATION/STATS

# List 1 contains the complete dataset for each network.
# df_DDI_tot contains all networks combined.
font_size=2

df_DDI_tot <- read_csv(file.path("data/5.Targets_NetworkDistance",
                                 "df_DDI_metrics_all.csv")) |> 
  mutate(int_sign_ebw=factor(int_sign_ebw,levels=c("Synergy","Additivity","Antagonism"))) 

# Minimum path ~ interaction type (ebw). 
my_comparisons <- list( c("Synergy", "Additivity"), c("Additivity", "Antagonism"), c("Synergy", "Antagonism") )
net <- df_DDI_tot$network |> unique() |> sort()

strains <- c("int_sign_ebw","int_sign_ecr","int_sign_seo","int_sign_stm","int_sign_pae","int_sign_pau")

for (j in seq_along(strains)){
  str_name <- strains[j]

plot_path <- list()
for (i in seq_along(net)){
  
  max_k <- df_DDI_tot |> filter(network==net[i]) |>
    filter(mean.path.length!="Inf") |> select(mean.path.length) |> max()
  
  min_k <- df_DDI_tot |> filter(network==net[i]) |>
    filter(mean.path.length!="Inf") |> select(mean.path.length) |> min()
  
  len_k <- max_k - min_k
  
  top <- max_k +len_k*0.3
  
  bot1 <- min_k - len_k*0.1
  
  bot2 <- min_k - len_k*0.2
  
  plot_path[[i]] <- df_DDI_tot |> filter(network==net[i]) |> # See also net[2]
    filter(mean.path.length!="Inf") |> 
    ggplot(aes(x=int_sign_ebw,y=mean.path.length,fill=int_sign_ebw))+
    geom_violin()+
    geom_boxplot(width=0.01,fill="white")+
    stat_compare_means(comparisons = my_comparisons,paired = FALSE)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = bot1) +
    stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
    theme_minimal()+
    scale_fill_brewer(palette="Set2")+
    xlab("Type of interaction (ebw)")+
    ylab("Minimum path length between connected protein targets")+
    theme(legend.position = "none")+
    ggtitle(net[i])+
    stat_n_text(y.pos = bot2)+
    ylim(bot2,top)
}

wrap_plots(plot_path[[1]],plot_path[[3]])


# Synergies are smaller than additive and antagonistic interactions in the PPI network.
# No signif difference between add/ant in LC

# Connectivity ~ interaction type (ebw).

plot_Kedge <- list()
for (i in seq_along(net)){
  max_k <- df_DDI_tot |> filter(network==net[i]) |>
    filter(mean.path.length!="Inf") |> select(mean.k.edge) |> max()
  
  min_k <- df_DDI_tot |> filter(network==net[i]) |>
    filter(mean.path.length!="Inf") |> select(mean.k.edge) |> min()
  
  len_k <- max_k - min_k
  
  top <- max_k +len_k*0.3
  
  bot1 <- min_k - len_k*0.1
  
  bot2 <- min_k - len_k*0.2
  
plot_Kedge[[i]] <- df_DDI_tot |> filter(network==net[i]) |>
    filter(mean.path.length!="Inf") |> 
    ggplot(aes(x=int_sign_ebw,y=mean.k.edge,fill=int_sign_ebw))+
    geom_violin()+
    geom_boxplot(width=0.01,fill="white")+
    stat_compare_means(comparisons = my_comparisons, paired=FALSE)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = bot1) +
    stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
    theme_minimal()+
    scale_fill_brewer(palette="Set2")+
    xlab("Type of interaction (ebw)")+
    ylab("K edge connectivity between connected protein targets")+
    theme(legend.position = "none")+
    ggtitle(net[i])+
    stat_n_text(y.pos = bot2)+
    ylim(bot2,top)
}

wrap_plots(plot_Kedge[[1]],plot_Kedge[[3]])
#Connectivity is lower in antagonism than additivities, and lower than synergies for co-functional.


# Node degree ~ interaction type (ebw). 
plot_deg <- list()
for (i in seq_along(net)){
  max_k <- df_DDI_tot |> filter(network==net[i]) |>
    filter(mean.path.length!="Inf") |> select(mean.mean.degree) |> max()
  
  min_k <- df_DDI_tot |> filter(network==net[i]) |>
    filter(mean.path.length!="Inf") |> select(mean.mean.degree) |> min()
  
  len_k <- max_k - min_k
  
  top <- max_k +len_k*0.3
  
  bot1 <- min_k - len_k*0.1
  
  bot2 <- min_k - len_k*0.2
  
  plot_deg[[i]] <- df_DDI_tot |> filter(network==net[i]) |>
    filter(mean.path.length!="Inf") |> 
    ggplot(aes(x=int_sign_ebw,y=mean.mean.degree,fill=int_sign_ebw))+
    geom_violin()+
    geom_boxplot(width=0.01,fill="white")+
    stat_compare_means(comparisons = my_comparisons, paired=FALSE)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = bot1) +
    stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
    theme_minimal()+
    scale_fill_brewer(palette="Set2")+
    xlab("Type of interaction (ebw)")+
    ylab("Average node degree between protein targets")+
    theme(legend.position = "none")+
    ggtitle(net[i])+
    stat_n_text(y.pos = bot2)+
    ylim(bot2,top)
}

wrap_plots(plot_deg[[1]],plot_deg[[3]])
# Significant differences in mean degree for co-functional (Eco-Cyc/GO-BP)


# Betweenness centrality ~ interaction type (ebw). 
plot_bet <- list()
for (i in seq_along(net)){
  max_k <- df_DDI_tot |> filter(network==net[i]) |>
    filter(mean.path.length!="Inf") |> select(mean.mean.bet) |> max()
  
  min_k <- df_DDI_tot |> filter(network==net[i]) |>
    filter(mean.path.length!="Inf") |> select(mean.mean.bet) |> min()
  
  len_k <- max_k - min_k
  
  top <- max_k +len_k*0.3
  
  bot1 <- min_k - len_k*0.1
  
  bot2 <- min_k - len_k*0.2
  
  plot_bet[[i]] <- df_DDI_tot |> filter(network==net[i]) |>
    filter(mean.path.length!="Inf") |> 
    ggplot(aes(x=int_sign_ebw,y=mean.mean.bet,fill=int_sign_ebw))+
    geom_violin()+
    geom_boxplot(width=0.01,fill="white")+
    stat_compare_means(comparisons = my_comparisons, paired=FALSE)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = bot1) +
    stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
    theme_minimal()+
    scale_fill_brewer(palette="Set2")+
    xlab("Type of interaction (ebw)")+
    ylab("Average node betweeness centrality between protein targets")+
    theme(legend.position = "none")+
    ggtitle(net[i])+
    stat_n_text(y.pos = bot2)+
    ylim(bot2,top)
}
wrap_plots(plot_bet[[1]],plot_bet[[3]])
# Significant differences in mean betweenness centrality for co-functional (Eco-Cyc/GO-BP)

# Eigenvector centrality ~ interaction type (ebw). 
plot_EVC <- list()
for (i in seq_along(net)){
  max_k <- df_DDI_tot |> filter(network==net[i]) |>
    filter(mean.path.length!="Inf") |> select(mean.mean.EVC) |> max()
  
  min_k <- df_DDI_tot |> filter(network==net[i]) |>
    filter(mean.path.length!="Inf") |> select(mean.mean.EVC) |> min()
  
  len_k <- max_k - min_k
  
  top <- max_k +len_k*0.3
  
  bot1 <- min_k - len_k*0.1
  
  bot2 <- min_k - len_k*0.2
  
  plot_EVC[[i]] <- df_DDI_tot |> filter(network==net[i]) |>
    filter(mean.path.length!="Inf") |> 
    ggplot(aes(x=int_sign_ebw,y=mean.mean.EVC,fill=int_sign_ebw))+
    geom_violin()+
    geom_boxplot(width=0.01,fill="white")+
    stat_compare_means(comparisons = my_comparisons, paired=FALSE)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = bot1) +
    stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
    theme_minimal()+
    scale_fill_brewer(palette="Set2")+
    xlab("Type of interaction (ebw)")+
    ylab("Average node eigenvector centrality between protein targets")+
    theme(legend.position = "none")+
    ggtitle(net[i])+
    stat_n_text(y.pos = bot2)+
    ylim(bot2,top)
}
wrap_plots(plot_EVC[[1]],plot_EVC[[3]])
}
# Significant differences in mean eigenvector centrality for co-functional (Eco-Cyc/GO-BP)


# Other visualizations:
# % synergies ~ network metric

# Total synergies ~ network metric





# % Connected ~ Type
# For each type of DDI, how many have its targets connected?
list_conn <- list()
for (i in 1:length(filenames)){
  list_conn[[i]] <- df_DDI_tot |>
    filter(network==net[i]) |> 
    select(drug_pair,connection_onoff,int_sign_ebw)  |> 
    group_by(int_sign_ebw,connection_onoff) |> 
    summarise(count=n()) |> 
    ggplot(aes(x = int_sign_ebw, fill = connection_onoff,y=count))+
    geom_col(position = "dodge")+
    theme_minimal()+
    scale_fill_brewer("Target connection",palette="Set2")+
    xlab("Type of interaction (ebw)")+
    ylab("Number of DDIs")+
    ggtitle(net[i])+
    theme(legend.position = "bottom")+
    geom_text(aes(label = paste0("n=",count)), vjust = 0, hjust = 0)+
    ylim(0,170)
}

wrap_plots(list_conn[[1]],list_conn[[3]])
#connected green, disconnected orange.
#Connected is more common than disconnected in all categories.
# All the antagonism found in metabolic network are disconnected.
# Even split in additivities (metabolic) for syn/ant

# Rate ~ Adj
#Adjacent proteins have a higher sigma rate than non-adjacent.

listAdj_met <- df_DDI_tot |>
  filter(network==net[5]) |> 
  mutate(max.adjacency=as.logical(max.adjacency)) |> 
  ggplot(aes(y=sigma.rate,x=as.factor(max.adjacency),
             fill=as.factor(max.adjacency)))+
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = -3,size=font_size) +
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
  theme_minimal()+
  theme(legend.position = "none")+
  ggtitle(net[5])+
  ylab(expression("Sigma rate ("~Bliss^2/MYA~")"))+
  xlab("Adjacent")+
  ylim(0,0.0090)+
  stat_n_text(y.pos = 0)

listAdj_PPI <- df_DDI_tot |>
  filter(network==net[3]) |> 
  mutate(max.adjacency=as.logical(max.adjacency)) |> 
  ggplot(aes(y=sigma.rate,x=as.factor(max.adjacency),
             fill=as.factor(max.adjacency)))+
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = -3,size=font_size) +
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
  theme_minimal()+
  theme(legend.position = "none")+
  ggtitle(net[3])+
  ylab(expression("Sigma rate ("~Bliss^2/MYA~")"))+
  xlab("Adjacent")+
  ylim(0,0.0090)+
  stat_n_text(y.pos = 0)


# Rate ~ Conn
#Connected proteins (EcoCyc, metabolic network) have a higher sigma rate than disconnected proteins.

my_comparisons <- list(c("FALSE","TRUE"))
listCon_met <- df_DDI_tot |>
  filter(network==net[5]) |> 
  mutate(connection_onoff=ifelse(connection_onoff=="connected","1","0")) |> 
  mutate(connection_onoff=as.logical(as.numeric(connection_onoff))) |> 
  ggplot(aes(y=sigma.rate,x=connection_onoff,
             fill=connection_onoff))+
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = -3,size=font_size) +
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
  theme_minimal()+
  theme(legend.position = "none")+
  ggtitle(net[5])+
  ylab(expression("Sigma rate ("~Bliss^2/MYA~")"))+
  xlab("Connected")+
  ylim(0,+0.0090)+
  stat_n_text(y.pos = 0)

listCon_PPI <- df_DDI_tot |>
  filter(network==net[3]) |> 
  mutate(connection_onoff=ifelse(connection_onoff=="connected","1","0")) |> 
  mutate(connection_onoff=as.logical(as.numeric(connection_onoff))) |> 
  ggplot(aes(y=sigma.rate,x=connection_onoff,
             fill=connection_onoff))+
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = -3,size=font_size) +
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
  theme_minimal()+
  theme(legend.position = "none")+
  ggtitle(net[3])+
  ylab(expression("Sigma rate ("~Bliss^2/MYA~")"))+
  xlab("Connected")+
  stat_n_text(y.pos = 0)

############################################
# For each cluster, calculate:

# number of unique DDIs per cluster
# df_DDI_tot  |> 
#   select(network,sigma.rate,clusters,drug_pair,connection_onoff,int_sign_ebw,sum_g) |> 
#   distinct() |> 
#   group_by(sigma.rate,network) |>
#   summarise(count = n_distinct(drug_pair)) |> 
#   ggplot(aes(x=sigma.rate,y=count,group=sigma.rate))+geom_point()+theme_minimal()+
#   facet_wrap(.~network)+
#   ylab("Number of DDIs")
# 
# # average sum_g
# df_DDI_tot |>
#   filter(network=="Small/medium-scale PPI") |> 
#   select(sigma.rate,clusters,drug_pair,connection_onoff,int_sign_ebw,sum_g,
#          sum_g)  |> 
#   distinct() |> 
#   ggplot(aes(x=sigma.rate,y=sum_g,group=sigma.rate))+geom_boxplot()+
#   geom_jitter()+theme_minimal()
# # clusters with higher rates are mainly a single int type (synergies or antagonism).


# Supplementary figure.
# SupFig3.Adjacency
sup_a <- listAdj_met
sup_b <- listAdj_PPI
# SupFig4.Connectedness
sup_c <- listCon_met
# SupFig5.PercentageConnected
sup_g <- list_conn[[3]]
sup_h <- list_conn[[11]]+
  theme(legend.position = "none")


# Evolutionary rates ~ int_type/distance/connectivity/node 

# Rates ~ Interaction Type (ebw) (PPI). 5C

df_target_tot <- left_join(full_df,df_target_tot,by="drugdrug")

c <- df_target_tot |>
  filter(network==net[1]|network==net[3]) |> 
  ggline(x = "int_sign_ebw.x", y = "sigma.rate.y", add = "mean_se",
         group = "network", col="network")+
  xlab("Type of interaction (ebw)")+
  ylab(expression("Sigma rate ("~Bliss^2/MYA~")"))+
  theme_minimal()+
  theme(legend.position = "none")
#Red:metabolic,blue:PPI
# Rates ~ path length. 5D
d <- df_target_tot |>
  filter(network==net[1]|network==net[3]) |> 
  mutate(path.length.corr=ifelse(is.infinite(path.length),"Not connected",
                                        ifelse(path.length==0,"0",
                                               ifelse(path.length==1,"1",
                                                      ifelse(path.length==2,"2",
                                                             ifelse(path.length==3,"3",
                                                                    ifelse(path.length==4,"4",
                                                                           ifelse(path.length==5,"5",
                                                                                  ifelse(path.length>=6,">6","Error"))))))))) |> 
  mutate(path.length.corr=factor(path.length.corr,levels=c("0","1","2","3","4","5",">6","Not connected","Error"))) |>
  ggline(x = "path.length.corr", y = "sigma.rate.x", add = c("mean_se"),
         group = "network", col="network",linetype = 1)+
  xlab("Path length")+
  ylab(expression("Sigma rate ("~Bliss^2/MYA~")"))+
  theme_minimal()+
  theme(legend.position = "bottom")


# Rates ~ Distance (all networks).
df_target_tot |>  mutate(path.length.corr=ifelse(path.length==0,"0",
                                                 ifelse(path.length==1,"1",
                                                        ifelse(path.length==2,"2",
                                                               ifelse(path.length==3,"3",
                                                                      ifelse(path.length==4,"4",
                                                                             ifelse(path.length==5,"5","6+"))))))) |>
  mutate(path.length.corr=factor(path.length.corr,levels=c("0","1","2","3","4","5","6+"))) |>
  ggline(x = "path.length.corr", y = "sigma.rate.y", add = "mean_se",
         group = "network", palette = "jco")+
  stat_compare_means(aes(group = network), label = "p.signif",
                     label.y = c(40, 40, 40))+
  facet_wrap(~network)+
  xlab("Minimum distance between targets")+
  ylab("Sigma rate")+
  theme_minimal()

# Rates ~ Interaction Type (ebw) (all networks).
df_target_tot |> ggline(x = "int_sign_ebw.x", y = "sigma.rate.y", add = "mean_se",
                     group = "network", palette = "jco")+
  stat_compare_means(aes(group = network), label = "p.signif",
                     label.y = c(40, 40, 40))+
  facet_wrap(~network)+
  xlab("Type of interaction")+
  ylab("Sigma rate")+
  theme_minimal()

# Rates ~ Interaction Type (all strains) (all networks).
df_target_tot |> ggline(x = "sum_g.y", y = "sigma.rate.y", add = "mean_se",
                     group = "network", palette = "jco")+
  stat_compare_means(aes(group = network), label = "p.signif",
                     label.y = c(40, 40, 40))+
  facet_wrap(~network)+
  xlab("Type of interaction")+
  ylab("Sigma rate")+
  theme_minimal()

# mean path ~ interaction (Syn,Ant,Add,Add-Ant,Ant-Syn,Add-Syn,Add-Ant-Syn)

# Mean.rate. 5A
a <- full_df |>
  select(type,sigma.rate)  |> 
  group_by(type) |> 
  summarise(Mean.rate=mean(sigma.rate),
            sem.rate=std.error(sigma.rate),
            count=n()) |> 
  ungroup() |> 
  ggplot(aes(x=type,y=Mean.rate))+
  geom_bar(position=position_dodge(), stat="identity",
           fill='#6191F2',width=0.5)+
  geom_errorbar(position=position_dodge(0.5),aes(ymin=Mean.rate-sem.rate, ymax=Mean.rate+sem.rate), width=.2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+
  xlab("Interaction type across all strains")+
  ylab(expression("Sigma rate ("~Bliss^2/MYA~")"))+
  theme(legend.position = "none")+
  geom_text(aes(label = paste0("n=",count)),vjust=-0.75)+
  ylim(0,0.008)
  


antagonism_score <- -1
synergy_score <- +1

b <- full_df |> 
  mutate(ebw_g=ifelse(int_sign_ebw=="Antagonism",antagonism_score,
                      ifelse(int_sign_ebw=="Additivity",0,
                             ifelse(int_sign_ebw=="Synergy",synergy_score,NA)))) |> 
  mutate(ecr_g=ifelse(int_sign_ecr=="Antagonism",antagonism_score,
                      ifelse(int_sign_ecr=="Additivity",0,
                             ifelse(int_sign_ecr=="Synergy",synergy_score,NA)))) |> 
  mutate(seo_g=ifelse(int_sign_seo=="Antagonism",antagonism_score,
                      ifelse(int_sign_seo=="Additivity",0,
                             ifelse(int_sign_seo=="Synergy",synergy_score,NA)))) |> 
  mutate(stm_g=ifelse(int_sign_stm=="Antagonism",antagonism_score,
                      ifelse(int_sign_stm=="Additivity",0,
                             ifelse(int_sign_stm=="Synergy",synergy_score,NA)))) |> 
  mutate(pae_g=ifelse(int_sign_pae=="Antagonism",antagonism_score,
                      ifelse(int_sign_pae=="Additivity",0,
                             ifelse(int_sign_pae=="Synergy",synergy_score,NA)))) |> 
  mutate(pau_g=ifelse(int_sign_pau=="Antagonism",antagonism_score,
                      ifelse(int_sign_pau=="Additivity",0,
                             ifelse(int_sign_pau=="Synergy",synergy_score,NA)))) |> 
  rowwise() |> 
  mutate(sum_g=ebw_g+ecr_g+seo_g+stm_g+pae_g+pau_g) |> 
  distinct()|> 
  group_by(sum_g) |> 
  summarise(Mean.rate=mean(sigma.rate),
            sem.rate=std.error(sigma.rate),
            count=n()) |> 
  ungroup() |> 
  ggplot(aes(x=sum_g,y=Mean.rate))+
  geom_bar(position=position_dodge(), stat="identity",
           fill='#6191F2',width=0.5)+
  geom_errorbar(position=position_dodge(0.5),aes(ymin=Mean.rate-sem.rate, ymax=Mean.rate+sem.rate), width=.2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+
  xlab("Score (More antagonistic <==> More synergistic)")+
  ylab(expression("Sigma rate ("~Bliss^2/MYA~")"))+
  theme(legend.position = "none")+
  geom_text(aes(label = paste0("n=",count)),vjust=-0.75)

# Mean.path. Sup_D
sup_d <- df_DDI_tot |>
  filter(network==net[1]|network==net[3]) |> 
  filter(!is.infinite(mean.path.length)) |> 
  select(network,sigma.rate,clusters,drug_pair,connection_onoff,type,mean.path.length)  |> 
  distinct() |> 
  group_by(type,network) |> 
  summarise(Mean.path=mean(mean.path.length),
            sem.path=std.error(mean.path.length),
            count=paste0("n=",n())) |> 
  ungroup() |> 
  ggplot(aes(x=type,y=Mean.path,fill=network))+
  geom_bar(position=position_dodge(), stat="identity",
           colour='black',width=0.5)+
  geom_errorbar(position=position_dodge(0.5),aes(ymin=Mean.path-sem.path, ymax=Mean.path+sem.path), width=.2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+
  xlab("Interaction type across all strains")+
  ylab("Minimum path distance of connected targets")+
  theme(legend.position = "none")+
  geom_text(aes(label = count,vjust=-0.75,hjust=1))+
  ylim(0,6)


#Mean.K.edge
sup_e <- df_DDI_tot |>
  filter(network==net[1]|network==net[3]) |> 
  select(network,sigma.rate,clusters,drug_pair,connection_onoff,type,mean.path.length,mean.k.edge)  |> 
  filter(!is.infinite(mean.path.length)) |> 
  distinct() |> 
  group_by(type,network) |> 
  summarise(Mean.K.edge=mean(mean.k.edge,na.rm=T),
            sem.K.edge=std.error(mean.k.edge),
            count=paste0("n=",n())) |> 
  ungroup() |> 
  ggplot(aes(x=type,y=Mean.K.edge,fill=network))+
  geom_bar(position=position_dodge(), stat="identity",
           colour='black',width=0.5)+
  geom_errorbar(position=position_dodge(0.5),aes(ymin=Mean.K.edge-sem.K.edge, ymax=Mean.K.edge+sem.K.edge), width=.2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+
  xlab("Interaction type across all strains")+
  ylab("Mean K-edge connectivity of connected targets")+
  theme(legend.position = "bottom")+
  geom_text(aes(label = count,vjust=-0.75,hjust=-0.25))+
  ylim(0,35)


# Mean Deg
sup_f <- df_DDI_tot |>
  filter(network==net[1]|network==net[3]) |> 
  select(network,sigma.rate,clusters,drug_pair,connection_onoff,type,mean.path.length,
         mean.mean.degree)  |> 
  distinct() |> 
  group_by(type,network) |> 
  summarise(Mean.mean_deg=mean(mean.mean.degree,na.rm=T),
            sem.mean_deg=std.error(mean.mean.degree),
            count=paste0("n=",n())) |> 
  ungroup() |> 
  ggplot(aes(x=type,y=Mean.mean_deg,fill=network))+
  geom_bar(position=position_dodge(), stat="identity",
           colour='black',width=0.5)+
  geom_errorbar(position=position_dodge(0.5),aes(ymin=Mean.mean_deg-sem.mean_deg, ymax=Mean.mean_deg+sem.mean_deg), width=.2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+
  xlab("Interaction type across all strains")+
  ylab("Mean node degree")+
  theme(legend.position = "none")+
  geom_text(aes(label = count,vjust=-0.75,hjust=-0.25))+
  ylim(0,45)


sup_i <- df_DDI_tot |>
  filter(network==net[1]|network==net[3]) |> 
  select(network,sigma.rate,clusters,drug_pair,connection_onoff,type,mean.path.length,
         mean.mean.degree,mean.k.edge) |> 
  distinct() |> 
  mutate(mean.k.group=ifelse(mean.k.edge<=5,"0-5",
                             ifelse(mean.k.edge<=10,"5-10",
                                    ifelse(mean.k.edge<=20,"10-20",
                                           ifelse(mean.k.edge>=20,">20","error"))))) |> 
  mutate(mean.k.group=factor(mean.k.group,levels=c("0-5","5-10","10-20",">20","error"))) |>
  ggline(x = "mean.k.group", y = "sigma.rate", add = "mean_se",
         group = "network", col="network")+
  xlab("K edge connectivity bettwen target proteins")+
  ylab(expression("Sigma rate ("~Bliss^2/MYA~")"))+
  theme_minimal()+
  theme(legend.position = "none")
  

sup_j <- df_DDI_tot |>
  filter(network==net[1]|network==net[3]) |> 
  select(network,sigma.rate,clusters,drug_pair,connection_onoff,type,mean.path.length,
         mean.mean.degree) |> 
  distinct() |> 
  mutate(mean.deg.group=ifelse(mean.mean.degree<=5,"0-5",
                             ifelse(mean.mean.degree<=10,"5-10",
                                    ifelse(mean.mean.degree<=20,"10-20",
                                           ifelse(mean.mean.degree>=20,">20","error"))))) |> 
  mutate(mean.deg.group=factor(mean.deg.group,levels=c("0-5","5-10","10-20",">20","error"))) |>
  ggline(x = "mean.deg.group", y = "sigma.rate", add = "mean_se",
         group = "network", col="network")+
  xlab("Average degree between target proteins")+
  ylab(expression("Sigma rate ("~Bliss^2/MYA~")"))+
  theme_minimal()+
  theme(legend.position = "none")

# Make figures:
# Fig 3
wrap_plots(plot_path_met,plot_path_PPI,
           plot_Kedge_met,plot_Kedge_PPI,
           plot_degree_met,plot_degree_PPI,byrow = FALSE)+ 
  plot_annotation(tag_levels = 'A',tag_suffix = '.')&
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 15, hjust = 0, vjust = 0))  #Fig3



# Fig 4
wrap_plots(a,b,c,d,nrow=2)+ 
  plot_annotation(tag_levels = 'A',tag_suffix = '.')&
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 15, hjust = 0, vjust = 0)) #Fig4


# 20 x 20. SupFig
wrap_plots(sup_a,sup_b,sup_c,sup_d,sup_e,sup_f,sup_g,
           sup_h,sup_i,sup_j,ncol=3)+ 
  plot_annotation(tag_levels = 'A',tag_suffix = '.')&
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 15, hjust = 0, vjust = 0)) 

df_DDI_tot |> 
  filter(network==net[1]|network==net[3]) |> 
  filter(Drug1 %in% c("A22","Novobiocin","Ciprofloxacin",
                      "Meropenem","Sulfamonomethoxine","Trimethoprim")) |>
  filter(Drug2 %in% c("A22","Novobiocin","Ciprofloxacin",
                      "Meropenem","Sulfamonomethoxine","Trimethoprim")) |> 
  select(network,drug_pair,clusters,sigma.rate,
         mean.path.length,mean.k.edge,mean.mean.degree,
         starts_with("int"),starts_with("SR")) |> 
  filter(drug_pair %in% c("Sulfamonomethoxine_Trimethoprim",
                          "A22_Novobiocin",
                          "Ciprofloxacin_Meropenem")) |> 
  arrange(drug_pair,network) |> 
  write.csv("data/5.Targets_NetworkDistance/sup_table1.csv",row.names = FALSE)



#Rates ~ interaction type and resistance type.
df_DDI_tot |> ggline(x = "int_sign_ebw", y = "sigma.rate", add = "mean_se",
                     group = "network", palette = "jco")+
  stat_compare_means(aes(group = network), label = "p.signif",
                     label.y = c(40, 40, 40))+
  facet_wrap(~network)+
  xlab("Type of interaction")+
  ylab("Sigma rate")+
  theme_minimal()


df1 <- df_DDI_tot |> 
  filter(network==net[1]|network==net[3]) |> 
  filter(drug_pair %in% c("A22_Novobiocin","Ciprofloxacin_Meropenem",
                          "Sulfamonomethoxine_Trimethoprim"))

