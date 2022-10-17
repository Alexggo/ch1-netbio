library(pacman)
p_load(igraph,tidyverse,BioNet,ape,geomorph,matrixStats,
       ComplexHeatmap,ggpubr,plotrix,tidymodels,patchwork,here)
library(EnvStats) # for adding sample size
source("bin/5.1.Networks.R")

# VISUALIZATION/STATS

# List 1 contains the complete dataset for each network.
# df_DDI_tot contains all networks combined.
#Execute after 5.1
font_size=2

df_DDI_tot <- read_csv(file.path("data/5.Targets_NetworkDistance","df_DDI_tot_metrics.csv")) |> 
  mutate(int_sign_ebw=factor(int_sign_ebw,levels=c("Synergy","Additivity","Antagonism"))) |> 
  arrange(network)

# Minimum path ~ interaction type (ebw). 
my_comparisons <- list( c("Synergy", "Additivity"), c("Additivity", "Antagonism"), c("Synergy", "Antagonism") )
net <- df_DDI_tot$network |> unique()

plot_path_met <- df_DDI_tot |> filter(network==net[3]) |> # See also net[5]
  filter(mean.path.length!="Inf") |> 
  ggplot(aes(x=int_sign_ebw,y=mean.path.length,fill=int_sign_ebw))+
  geom_violin()+
  geom_boxplot(width=0.01,fill="white")+
  stat_compare_means(comparisons = my_comparisons,paired = FALSE)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.5) +
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
  theme_minimal()+
  scale_fill_brewer(palette="Set2")+
  xlab("Type of interaction (ebw)")+
  ylab("Minimum path length between connected protein targets")+
  theme(legend.position = "none")+
  ggtitle(net[3])+
  ylim(c(0,7.5))+
  stat_n_text(y.pos = 0)

plot_path_PPI<- df_DDI_tot |> filter(network==net[11]) |> 
  filter(mean.path.length!="Inf") |> 
  ggplot(aes(x=int_sign_ebw,y=mean.path.length,fill=int_sign_ebw))+
  geom_violin()+
  geom_boxplot(width=0.01,fill="white")+
  stat_compare_means(comparisons = my_comparisons,paired = FALSE)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.5) +
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
  theme_minimal()+
  scale_fill_brewer(palette="Set2")+
  xlab("Type of interaction (ebw)")+
  ylab("Minimum path length between connected protein targets")+
  theme(legend.position = "none")+
  ggtitle(net[11])+
  ylim(c(0,7.5))+
  stat_n_text(y.pos = 0)

wrap_plots(plot_path_met,plot_path_PPI)


# Synergies are smaller than additive and antagonistic interactions in the PPI network.
# No signif difference between add/ant in LC

# Connectivity ~ interaction type (ebw).
my_comparisons <- list( c("Synergy", "Additivity"), c("Additivity", "Antagonism"), c("Synergy", "Antagonism") )
net <- df_DDI_tot$network |> unique()

plot_Kedge_met <- df_DDI_tot |> filter(network==net[3]) |>
  filter(mean.path.length!="Inf") |> 
  ggplot(aes(x=int_sign_ebw,y=mean.k.edge,fill=int_sign_ebw))+
  geom_violin()+
  geom_boxplot(width=0.01,fill="white")+
  stat_compare_means(comparisons = my_comparisons, paired=FALSE)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = -3) +
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
  theme_minimal()+
  scale_fill_brewer(palette="Set2")+
  xlab("Type of interaction (ebw)")+
  ylab("K edge connectivity between connected protein targets")+
  theme(legend.position = "none")+
  ggtitle(net[3])+
  stat_n_text(y.pos = -8)+
  ylim(-10,65)

plot_Kedge_PPI <- df_DDI_tot |> filter(network==net[11]) |>
  filter(mean.path.length!="Inf") |> 
  ggplot(aes(x=int_sign_ebw,y=mean.k.edge,fill=int_sign_ebw))+
  geom_violin()+
  geom_boxplot(width=0.01,fill="white")+
  stat_compare_means(comparisons = my_comparisons, paired=FALSE)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0) +
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
  theme_minimal()+
  scale_fill_brewer(palette="Set2")+
  xlab("Type of interaction (ebw)")+
  ylab("K edge connectivity between connected protein targets")+
  theme(legend.position = "none")+
  ggtitle(net[11])+
  stat_n_text(y.pos = -0.4)+
  ylim(-.5,5)

wrap_plots(plot_Kedge_met,plot_Kedge_PPI)

#Connectivity is lower in antagonism than additivities, and lower than synergies for co-functional.


# Node degree ~ interaction type (ebw). 
my_comparisons <- list( c("Synergy", "Additivity"), c("Additivity", "Antagonism"), c("Synergy", "Antagonism") )
net <- df_DDI_tot$network |> unique()

plot_degree_met <- df_DDI_tot |> filter(network==net[3]) |> 
  ggplot(aes(x = int_sign_ebw, y = mean.mean.degree,
             fill = int_sign_ebw))+ 
  geom_violin()+
  geom_boxplot(width=0.01,fill="white")+
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0) +
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
  theme_minimal()+
  scale_fill_brewer(palette="Set2")+
  xlab("Type of interaction (ebw)")+
  ylab("Average node degree between protein targets")+
  theme(legend.position = "none")+
  ggtitle(net[3])+
  stat_n_text(y.pos = -3)

plot_degree_PPI <- df_DDI_tot |> filter(network==net[11]) |> 
  ggplot(aes(x = int_sign_ebw, y = mean.mean.degree,
             fill = int_sign_ebw))+ 
  geom_violin()+
  geom_boxplot(width=0.01,fill="white")+
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = -0) +
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
  theme_minimal()+
  scale_fill_brewer(palette="Set2")+
  xlab("Type of interaction (ebw)")+
  ylab("Average node degree between protein targets")+
  theme(legend.position = "none")+
  ggtitle(net[11])+
  stat_n_text(y.pos = -0.5)+
  ylim(-0.5,8)

wrap_plots(plot_degree_met,plot_degree_PPI)
# Significant differences in mean degree for co-functional (Eco-Cyc/GO-BP)

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
    ylab("Number of connections")+
    ggtitle(net[i])+
    theme(legend.position = "bottom")+
    geom_text(aes(label = paste0("n=",count)), vjust = 0, hjust = 0)
}

wrap_plots(list_conn[[3]],list_conn[[11]])
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
  ylim(0,0.0080)+
  stat_n_text(y.pos = 0)

listAdj_PPI <- df_DDI_tot |>
  filter(network==net[11]) |> 
  mutate(max.adjacency=as.logical(max.adjacency)) |> 
  ggplot(aes(y=sigma.rate,x=as.factor(max.adjacency),
             fill=as.factor(max.adjacency)))+
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = -3,size=font_size) +
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
  theme_minimal()+
  theme(legend.position = "none")+
  ggtitle(net[11])+
  ylab(expression("Sigma rate ("~Bliss^2/MYA~")"))+
  xlab("Adjacent")+
  ylim(0,0.0080)+
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
  ylim(0,+0.0080)+
  stat_n_text(y.pos = 0)

listCon_PPI <- df_DDI_tot |>
  filter(network==net[11]) |> 
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
  ggtitle(net[11])+
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
c <- df_target_tot |>
  filter(network==net[3]|network==net[11]) |> 
  ggline(x = "int_sign_ebw", y = "sigma.rate", add = "mean_se",
         group = "network", col="network")+
  xlab("Type of interaction (ebw)")+
  ylab(expression("Sigma rate ("~Bliss^2/MYA~")"))+
  theme_minimal()+
  theme(legend.position = "none")
#Red:metabolic,blue:PPI
# Rates ~ path length. 5D
d <- df_target_tot |>
  filter(network==net[3]|network==net[11]) |> 
  mutate(path.length.corr=ifelse(is.infinite(path.length),"Not connected",
                                        ifelse(path.length==0,"0",
                                               ifelse(path.length==1,"1",
                                                      ifelse(path.length==2,"2",
                                                             ifelse(path.length==3,"3",
                                                                    ifelse(path.length==4,"4",
                                                                           ifelse(path.length==5,"5",
                                                                                  ifelse(path.length>=6,">6","Error"))))))))) |> 
  mutate(path.length.corr=factor(path.length.corr,levels=c("0","1","2","3","4","5",">6","Not connected","Error"))) |>
  ggline(x = "path.length.corr", y = "sigma.rate", add = c("mean_se"),
         group = "network", col="network",linetype = 1)+
  xlab("Path length")+
  ylab(expression("Sigma rate ("~Bliss^2/MYA~")"))+
  theme_minimal()+
  theme(legend.position = "bottom")


r <- c()
for (i in 1:dim(df_DDI_tot)[1]){
  x <- df_DDI_tot[i,25:30] |> t() |> as.vector()
  x <- x |> unique() |> sort()
  r[i] <- paste0(x,collapse = "-")
}
df_DDI_tot$type <- r

# Rates ~ Distance (all networks).
df_target_tot |>  mutate(path.length.corr=ifelse(path.length==0,"0",
                                                 ifelse(path.length==1,"1",
                                                        ifelse(path.length==2,"2",
                                                               ifelse(path.length==3,"3",
                                                                      ifelse(path.length==4,"4",
                                                                             ifelse(path.length==5,"5","6+"))))))) |>
  mutate(path.length.corr=factor(path.length.corr,levels=c("0","1","2","3","4","5","6+"))) |>
  ggline(x = "path.length.corr", y = "sigma.rate", add = "mean_se",
         group = "network", palette = "jco")+
  stat_compare_means(aes(group = network), label = "p.signif",
                     label.y = c(40, 40, 40))+
  facet_wrap(~network)+
  xlab("Minimum distance between targets")+
  ylab("Sigma rate")+
  theme_minimal()

# Rates ~ Interaction Type (ebw) (all networks).
df_DDI_tot |> ggline(x = "int_sign_ebw", y = "sigma.rate", add = "mean_se",
                     group = "network", palette = "jco")+
  stat_compare_means(aes(group = network), label = "p.signif",
                     label.y = c(40, 40, 40))+
  facet_wrap(~network)+
  xlab("Type of interaction")+
  ylab("Sigma rate")+
  theme_minimal()

# Rates ~ Interaction Type (all strains) (all networks).
df_DDI_tot |> ggline(x = "sum_g", y = "sigma.rate", add = "mean_se",
                     group = "network", palette = "jco")+
  stat_compare_means(aes(group = network), label = "p.signif",
                     label.y = c(40, 40, 40))+
  facet_wrap(~network)+
  xlab("Type of interaction")+
  ylab("Sigma rate")+
  theme_minimal()

# mean path ~ interaction (Syn,Ant,Add,Add-Ant,Ant-Syn,Add-Syn,Add-Ant-Syn)

r <- c()
for (i in 1:dim(rates)[1]){
  x <- rates[i,16:21] |> t() |> as.vector()
  x <- x |> unique() |> sort()
  r[i] <- paste0(x,collapse = "-")
}
rates$type <- r

# Mean.rate. 5A
a <- rates |>
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

b <- rates |> 
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
  filter(network==net[3]|network==net[11]) |> 
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
  filter(network==net[3]|network==net[11]) |> 
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
  geom_text(aes(label = count,vjust=-0.75,hjust=-0.5))


# Mean Deg
sup_f <- df_DDI_tot |>
  filter(network==net[3]|network==net[11]) |> 
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
  geom_text(aes(label = count,vjust=-0.75,hjust=-0.5))


sup_i <- df_DDI_tot |>
  filter(network==net[3]|network==net[11]) |> 
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
  filter(network==net[3]|network==net[11]) |> 
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

df_DDI_tot |> 
  filter(clusters=="13") |> 
  filter(network==net[10]) |> 
  View()

# Make figures:
# Fig 3
wrap_plots(plot_path_met,plot_path_PPI,
           plot_Kedge_met,plot_Kedge_PPI,
           plot_degree_met,plot_degree_PPI,byrow = FALSE)+ 
  plot_annotation(tag_levels = 'A',tag_suffix = '.')&
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 15, hjust = 0, vjust = 0))  #Fig3



# Fig 5
wrap_plots(a,b,c,d,nrow=2)+ 
  plot_annotation(tag_levels = 'A',tag_suffix = '.')&
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 15, hjust = 0, vjust = 0)) #Fig5


# 20 x 20. SupFig
wrap_plots(sup_a,sup_b,sup_c,sup_d,sup_e,sup_f,sup_g,
           sup_h,sup_i,sup_j,ncol=3)+ 
  plot_annotation(tag_levels = 'A',tag_suffix = '.')&
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 15, hjust = 0, vjust = 0)) 

df_DDI_tot |> 
  filter(network==net[3]|network==net[11]) |> 
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
