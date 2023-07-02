# "5.2.Visualizations"
library(pacman)
p_load(igraph,tidyverse,BioNet,ape,geomorph,matrixStats,
       ggpubr,plotrix,tidymodels,patchwork,here,EnvStats,report)

set_name <- "sen2in1" #allddi, sen2in1
set_col <- "all" #all. #Run with the values for all in both allddi and sen2in
# VISUALIZATION/STATS

# List 1 contains the complete dataset for each network.
# df_DDI_tot contains all networks combined.
font_size=2

df_target_tot <- read_csv(paste0("results/","df_target_metrics",".csv"))

df_DDI_tot <- read_csv(paste0("results/","df_DDI_metrics",".csv")) |> 
  mutate(int_sign_ebw=factor(int_sign_ebw,levels=c("Synergy","Additivity","Antagonism","Non-targets"))) |> 
  mutate(int_sign_ecr=factor(int_sign_ecr,levels=c("Synergy","Additivity","Antagonism","Non-targets"))) |> 
  mutate(int_sign_stm=factor(int_sign_stm,levels=c("Synergy","Additivity","Antagonism","Non-targets"))) |> 
  mutate(int_sign_seo=factor(int_sign_seo,levels=c("Synergy","Additivity","Antagonism","Non-targets"))) |> 
  mutate(int_sign_pae=factor(int_sign_pae,levels=c("Synergy","Additivity","Antagonism","Non-targets"))) |> 
  mutate(int_sign_pau=factor(int_sign_pau,levels=c("Synergy","Additivity","Antagonism","Non-targets"))) |> 
  mutate(connection_onoff=ifelse(mean.k.edge==0,"FALSE","TRUE")) |> 
  mutate(connection_onoff=factor(connection_onoff,levels = c("TRUE","FALSE")))

# Figures Network metric \~ interaction score.

# Minimum path ~ interaction type (ebw). 
my_comparisons <- list( c("Synergy", "Additivity"), c("Additivity", "Antagonism"), c("Synergy", "Antagonism"),c("Synergy", "Non-targets"), c("Additivity", "Non-targets"), c("Antagonism", "Non-targets"))
net <- df_DDI_tot$network |> unique() |> sort()
strains <- c("int_sign_ebw","int_sign_ecr","int_sign_seo","int_sign_stm","int_sign_pae","int_sign_pau")
short_str <-   str_split(strains,"_") |> map_chr(`[`,3)

y_axis <- c("Minimum path length between connected protein targets",
            "K edge connectivity between connected protein targets",
            "Average node degree",
            "Average betweeness centrality",
            "Average eigenvector centrality",
            "Average sigma rate","Average sigma rate (set2in1)")
y_var <- c("mean.path.length","mean.k.edge","mean.mean.degree",
           "mean.mean.bet","mean.mean.EVC","mean.sigma_all","mean.sigma_set")


# set_sel is NA for the sample non-target nodes.
# set_sel is TRUE if the DDI is in the SET, and FALSE if not.
# set_sel = TRUE, or set_sel = c(TRUE,FALSE)
plot_net_int <- function(net_sel,var_sel,str_name,set_sel){
  max_k <- df_DDI_tot |> 
    filter(network==net[net_sel]) |> 
    filter(mean.path.length!="Inf") |> 
    filter(is.na(SET)|SET %in% set_sel) |> 
    select(y_var[var_sel]) |> 
    max()
  
  min_k <- df_DDI_tot |> 
    filter(network==net[net_sel]) |> 
    filter(is.na(SET)|SET %in% set_sel) |> 
    filter(mean.path.length!="Inf") |> 
    select(y_var[var_sel]) |> 
    min()
  
  len_k <- max_k - min_k
  top <- max_k +len_k*0.7
  bot1 <- min_k - len_k*0.1
  bot2 <- min_k - len_k*0.2
  net_name <- net[net_sel]
  
  df1 <- df_DDI_tot |> 
    filter(network==net[net_sel]) |> 
    filter(is.na(SET)|SET %in% set_sel) |>
    select(strains[str_name],y_var[var_sel],network)  |> 
    filter(y_var[var_sel]!="Inf")
  
  colnames(df1) <- c("int_type","variable_name","network")
  
  net_name <- df1$network |> unique()
  
  pl1 <- df1 |> 
    ggplot(aes(x=int_type,y=variable_name,
               fill=int_type))+
    xlab(paste0("Type of interaction (",short_str[str_name],")"))+
    ylab(y_axis[var_sel])+
    geom_violin()+
    geom_boxplot(width=0.01,fill="white")+
    stat_compare_means(comparisons = my_comparisons,paired = FALSE)+ # Add p-value
    stat_compare_means(label.y = bot1) +
    stat_summary(fun=mean, geom="point", shape=23, size=4, col="black",fill="red")+
    theme_minimal()+
    theme(legend.position = "none")+
    ggtitle(label  = element_text(net_name))+
    stat_n_text(y.pos = bot2)+
    ylim(bot2,top)+
    scale_fill_brewer(palette = "Set2")+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0))
  
  plot(pl1)
}


# For all DDIs.
# Results per strain:
str_fig <- 1
fig3_str <- wrap_plots(plot_net_int(1,1,str_fig,c(T,F)),plot_net_int(3,1,str_fig,c(T,F)),
                       plot_net_int(1,2,str_fig,c(T,F)),plot_net_int(3,2,str_fig,c(T,F)),
                       plot_net_int(1,3,str_fig,c(T,F)),plot_net_int(3,3,str_fig,c(T,F)),
                       plot_net_int(1,4,str_fig,c(T,F)),plot_net_int(3,4,str_fig,c(T,F)),
                       plot_net_int(1,5,str_fig,c(T,F)),plot_net_int(3,5,str_fig,c(T,F)),
                       plot_net_int(1,6,str_fig,c(T,F)),plot_net_int(3,6,str_fig,c(T,F)),
                       plot_net_int(1,7,str_fig,c(T,F)),plot_net_int(3,7,str_fig,c(T,F)),
                       nrow = 2, byrow = FALSE)

var_fig <- 1
fig3_var <- wrap_plots(plot_net_int(1,var_fig,1,c(T,F)),plot_net_int(3,var_fig,1,c(T,F)),
                       plot_net_int(1,var_fig,2,c(T,F)),plot_net_int(3,var_fig,2,c(T,F)),
                       plot_net_int(1,var_fig,3,c(T,F)),plot_net_int(3,var_fig,3,c(T,F)),
                       plot_net_int(1,var_fig,4,c(T,F)),plot_net_int(3,var_fig,4,c(T,F)),
                       plot_net_int(1,var_fig,5,c(T,F)),plot_net_int(3,var_fig,5,c(T,F)),
                       plot_net_int(1,var_fig,6,c(T,F)),plot_net_int(3,var_fig,6,c(T,F)),
                       nrow = 2, byrow = FALSE)

fig3 <- wrap_plots(plot_net_int(1,1,str_fig,c(T,F)),plot_net_int(3,1,str_fig,c(T,F)),
                   plot_net_int(1,2,str_fig,c(T,F)),plot_net_int(3,2,str_fig,c(T,F)),
                   plot_net_int(1,3,str_fig,c(T,F)),plot_net_int(3,3,str_fig,c(T,F)),
                   plot_net_int(1,4,str_fig,c(T,F)),plot_net_int(3,4,str_fig,c(T,F)),plot_net_int(1,5,str_fig,c(T,F)),plot_net_int(3,5,str_fig,c(T,F)),
                   ncol = 4)+ 
  plot_annotation(tag_levels = 'A',tag_suffix = '.')&
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 15, hjust = 0, vjust = 0))+
  theme(axis.title = element_text(size = 12),  # Set axis title size
        axis.text = element_text(size = 11))

# Fig 3 includes a category for non-targets (NAs in column set)
# P-values with non-targets vary due to differences between samples in every iteration.
# The other tests are consistent.
ggsave(
  filename = paste0("results/allddi/fig3_allddi.pdf"),
  plot = fig3,
  width = 15,
  height = 15
)

ggsave(
  filename = paste0("results/allddi/fig3_allddi.png"),
  plot = fig3,
  width = 15,
  height = 15,dpi=500
)

# Set takes only the drugs that are part of the sensitive group for comparison
fig3_set <- wrap_plots(plot_net_int(1,1,str_fig,c(T)),plot_net_int(3,1,str_fig,c(T)),
                       plot_net_int(1,2,str_fig,c(T)),plot_net_int(3,2,str_fig,c(T)),
                       plot_net_int(1,3,str_fig,c(T)),plot_net_int(3,3,str_fig,c(T)),
                       plot_net_int(1,4,str_fig,c(T)),plot_net_int(3,4,str_fig,c(T)),plot_net_int(1,5,str_fig,c(T)),plot_net_int(3,5,str_fig,c(T)),
                       ncol = 4)+ 
  plot_annotation(tag_levels = 'A',tag_suffix = '.')&
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 15, hjust = 0, vjust = 0))+
  theme(axis.title = element_text(size = 12),  # Set axis title size
        axis.text = element_text(size = 11))

ggsave(
  filename = paste0("results/sen2in1/fig3_sen2in1.pdf"),
  plot = fig3_set,
  width = 15,
  height = 15
)

ggsave(
  filename = paste0("results/sen2in1/fig3_sen2in1.png"),
  plot = fig3_set,
  width = 15,
  height = 15,dpi=500
)

# Relative number of DDI for each network metric
net <- df_DDI_tot$network |> unique() |> sort()
strains <- c("int_sign_ebw","int_sign_ecr","int_sign_seo","int_sign_stm","int_sign_pae","int_sign_pau")
short_str <-   str_split(strains,"_") |> map_chr(`[`,3)

x_axis <- c("Minimum path length between connected protein targets",
            "K edge connectivity between connected protein targets",
            "Average node degree",
            "Average betweeness centrality",
            "Average eigenvector centrality",
            "Average sigma rate","Average sigma rate (set2in1)")
x_var <- c("mean.path.length","mean.k.edge","mean.mean.degree",
           "mean.mean.bet","mean.mean.EVC","mean.sigma_all","mean.sigma_set")


plot_net_num <- function(net_sel,var_sel,str_name){
  max_k <- df_DDI_tot |> 
    filter(network==net[net_sel]) |> 
    filter(mean.path.length!="Inf") |> 
    select(y_var[var_sel]) |> 
    max()
  
  min_k <- df_DDI_tot |> 
    filter(network==net[net_sel]) |>
    filter(mean.path.length!="Inf") |> 
    select(y_var[var_sel]) |> 
    min()
  
  len_k <- max_k - min_k
  top <- max_k +len_k*0.5
  bot1 <- min_k - len_k*0.1
  bot2 <- min_k - len_k*0.2
  
  df1 <- df_DDI_tot |> filter(network==net[net_sel]) |> 
    filter(mean.path.length!="Inf") |>
    select(strains[str_name],y_var[var_sel],network)  
  
  colnames(df1)[1] <- "int_type"
  colnames(df1)[2] <- "variable_name"
  
  net_name <- df1$network |> unique()
  
  br_df1 <- df1 |> 
    mutate(variable_name =cut(variable_name,breaks=5)) |> 
    arrange(variable_name) |> 
    group_by(variable_name,network,int_type) |> 
    summarise(total=n()) |> 
    pivot_wider(names_from = int_type,values_from = total) |> 
    mutate(Synergy=ifelse(is.na(Synergy),0,Synergy))|> 
    mutate(Antagonism=ifelse(is.na(Antagonism),0,Antagonism))|> 
    mutate(Additivity=ifelse(is.na(Additivity),0,Additivity)) |> 
    mutate(Additivity=ifelse(is.na(`Non-targets`),0,`Non-targets`)) |> 
    rowwise() |> 
    mutate(total=Synergy+Antagonism+Additivity) |> 
    mutate(Perc.Synergy=round(100*Synergy/total,2),
           Perc.Antagonism=round(100*Antagonism/total,2),
           Perc.Additivity=round(100*Additivity/total,2)) |> 
    pivot_longer(names_to = "int_type",values_to = "percentage",8:10) |> 
    mutate(int_type=ifelse(int_type=="Perc.Synergy","Synergy",
                           ifelse(int_type=="Perc.Antagonism","Antagonism",
                                  ifelse(int_type=="Perc.Additivity","Additivity","NA")))) |>
    mutate(variable_name=factor(variable_name))
  
  pl1 <- br_df1 |> 
    ggplot(aes(x=variable_name,y=percentage,col=int_type,
               group=int_type))+
    xlab(y_axis[var_sel])+
    ylab("Percentage of DDI")+
    geom_point()+
    geom_line()+
    theme_minimal()+
    theme(legend.position = "bottom")+
    ggtitle(label  = element_text(net_name))+ 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  plot(pl1)
}

str_fig <- 1
fig2_str <- wrap_plots(plot_net_num(1,1,str_fig),plot_net_num(3,1,str_fig),
                       plot_net_num(1,2,str_fig),plot_net_num(3,2,str_fig),
                       plot_net_num(1,3,str_fig),plot_net_num(3,3,str_fig),
                       plot_net_num(1,4,str_fig),plot_net_num(3,4,str_fig),
                       plot_net_num(1,5,str_fig),plot_net_num(3,5,str_fig),
                       plot_net_num(1,6,str_fig),plot_net_num(3,6,str_fig),
                       plot_net_num(1,7,str_fig),plot_net_num(3,7,str_fig),
                       nrow = 2, byrow = FALSE)

var_fig <- 1
fig2_var <- wrap_plots(plot_net_num(1,var_fig,1),plot_net_num(3,var_fig,1),
                       plot_net_num(1,var_fig,2),plot_net_num(3,var_fig,2),
                       plot_net_num(1,var_fig,3),plot_net_num(3,var_fig,3),
                       plot_net_num(1,var_fig,4),plot_net_num(3,var_fig,4),
                       plot_net_num(1,var_fig,5),plot_net_num(3,var_fig,5),
                       plot_net_num(1,var_fig,6),plot_net_num(3,var_fig,6),
                       nrow = 2, byrow = FALSE)

# Comparisons and trends of network metrics. Fig 4 and Sup. Fig.

## Path length: Synergies are closer in the met and PPI network
#than additive and antagonistic interactions. 

## Connectivity: Synergies are more connected in the met network
#than additive and antagonistic interactions. 

## Degree: Synergies have a higher number of connections in the 
#met network than additivities/antagonisms. 

## Betweeness centrality: Synergies are more central in the met 
#network than antagonism and additivities. 

## Eigenvector centrality: Synergies are more central in the met
#network than antagonism and additivities. 

## Sigma squared per DDI: Synergies have greater sigma in the met
#and PPI networks squared values than antagonisms/additivities.


# % Connected ~ Type
# For each type of DDI, how many have its targets connected?
net <- df_DDI_tot$network |> unique() |> sort()
strains <- c("int_sign_ebw","int_sign_ecr","int_sign_seo","int_sign_stm","int_sign_pae","int_sign_pau")
short_str <-   str_split(strains,"_") |> map_chr(`[`,3)

conn_type <- function(net_sel,str_name){
  str_var <- strains[str_name]
  df2 <- df_DDI_tot |>
    filter(network==net[net_sel]) |> 
    select(DRUG_ID,connection_onoff,str_var)
  
  colnames(df2)[3] <- "str_var"
  pl2 <- df2 |>  group_by(str_var,connection_onoff) |> 
    summarise(count=n()) |> 
    ggplot(aes(x = str_var, fill = connection_onoff,y=count))+
    geom_col(position = "dodge")+
    theme_minimal()+
    scale_fill_brewer("Target connection",palette="Set2")+
    xlab(paste0("Type of interaction (",short_str[str_name],")"))+
    ylab("Number of DDIs")+
    ggtitle(net[net_sel])+
    theme(legend.position = "bottom")+
    geom_text(aes(label = paste0("n=",count)), vjust = 0, hjust = 0)
  
  plot(pl2)
}

wrap_plots(conn_type(1,1),conn_type(3,1))
wrap_plots(conn_type(1,2),conn_type(3,2))
wrap_plots(conn_type(1,3),conn_type(3,3))
wrap_plots(conn_type(1,4),conn_type(3,4))
wrap_plots(conn_type(1,5),conn_type(3,5))
wrap_plots(conn_type(1,6),conn_type(3,6))

#connected green, disconnected orange.
#Connected is more common than disconnected in all categories.
# All the antagonism found in metabolic network are disconnected.
# Even split in additivities (metabolic) for syn/ant

# Rate ~ Adj
my_comparisons <- list( c("FALSE", "TRUE"))
y_var <- c("max.adjacency","connection_onoff")
rate_adjcon <- function(net_sel,var_sel){
  var_1 <- y_var[var_sel]
  
  df3 <- df_DDI_tot |>
    filter(network==net[net_sel]) |> 
    filter(!is.na(SET)) |> 
    select(sigma.rate_all,max.adjacency,connection_onoff,network) |> 
    mutate(max.adjacency=as.logical(max.adjacency)) |> 
    mutate(max.adjacency=as.character(as.logical(max.adjacency))) |> 
    mutate(connection_onoff=as.character(connection_onoff)) |> 
    select(sigma.rate_all,var_1,network)
  
  colnames(df3)[2] <- "var_name"
  
  max_k <- df3  |> 
    select(sigma.rate_all) |> 
    max()
  
  min_k <- df3  |> 
    select(sigma.rate_all) |> 
    min()
  
  len_k <- max_k - min_k
  top <- max_k +len_k*0.5
  bot1 <- min_k - len_k*0.1
  bot2 <- min_k - len_k*0.2
  
  pl3 <- df3 |> 
    ggplot(aes(y=sigma.rate_all,x=var_name,
               fill=var_name))+
    geom_boxplot()+
    stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = bot1) +
    stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red")+
    theme_minimal()+
    theme(legend.position = "none")+
    ggtitle(net[net_sel])+
    ylab(expression("Sigma rate ("~Bliss^2/MYA~")"))+
    xlab(var_1)+
    stat_n_text(y.pos = bot2)+
    ylim(bot2,top)
  
  plot(pl3)
}

#Adjacent proteins have a higher sigma rate than non-adjacent.
rate_adjcon(1,1)
rate_adjcon(2,1)
rate_adjcon(3,1)

# Rate ~ Conn
#Connected proteins (EcoCyc, metabolic network) have a higher sigma rate than disconnected proteins.
rate_adjcon(1,2)
rate_adjcon(2,2)
rate_adjcon(3,2)

# Supplementary figure.
# SupFig.Adjacency
sup_a <- rate_adjcon(2,1)
sup_b <- rate_adjcon(3,1)
# SupFigc.Connectedness
sup_c <- rate_adjcon(2,2)
# SupFig i.Percentage Connected
sup_i <- conn_type(1,1)+
  theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))
sup_j <- conn_type(3,1)+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))

# Evolutionary rates ~ int_type/distance/connectivity/node 
# Fixed below this
df_target_DDI <- left_join(df_target_tot,df_DDI_tot,by="DRUG_ID",
                           suffix = c(".target", ".DDI")) |> 
  mutate(network=network.target) |> 
  filter(network==net[1]|network==net[3])

df <- read_csv(paste0("results/",set_name,"/","DDI_table_rates_",set_name,".csv"))
# mean path ~ interaction (Syn,Ant,Add,Add-Ant,Ant-Syn,Add-Syn,Add-Ant-Syn)

# Mean.rate. Fig4
fig.4a.df <- df |>
  select(type,sigma.rate) |> 
  filter(!is.na(type)) 

type_ddi <- fig.4a.df$type |> unique()
type_ddi1 <- fig.4a.df |> filter(type==type_ddi[1]) |> select(sigma.rate) |> pull()
type_ddi2 <- fig.4a.df |> filter(type==type_ddi[2]) |> select(sigma.rate) |> pull()
type_ddi3 <- fig.4a.df |> filter(type==type_ddi[3]) |> select(sigma.rate) |> pull()
type_ddi4 <- fig.4a.df |> filter(type==type_ddi[4]) |> select(sigma.rate) |> pull()
type_ddi5 <- fig.4a.df |> filter(type==type_ddi[5]) |> select(sigma.rate) |> pull()
type_ddi6 <- fig.4a.df |> filter(type==type_ddi[6]) |> select(sigma.rate) |> pull()

t_1_2 <- t.test(type_ddi1,type_ddi2)
t_1_2 |> glance()
t_1_2 |> report() #p < 0.002
t_1_3 <- t.test(type_ddi1,type_ddi3) 
t_1_3 |> glance()
t_1_3 |> report() #p < .001
t_1_4 <- t.test(type_ddi1,type_ddi4)
t_1_4 |> glance()
t_1_4 |> report() #p < .001
t_1_5 <- t.test(type_ddi1,type_ddi5)
t_1_5 |> glance()
t_1_5 |> report() #p < .001
#t_1_6 <- t.test(type_ddi1,type_ddi6) 
#t_1_6 |> glance()
#t_1_6 |> report() # Not enough observations
t_2_3 <- t.test(type_ddi2,type_ddi3)
t_2_3 |> glance()
t_2_3 |> report() #p < .001
t_2_4 <- t.test(type_ddi2,type_ddi4)
t_2_4 |> glance()
t_2_4 |> report() #p < 0.526. "Additivity-Antagonism-Synergy" vs "Additivity-Antagonism"
t_2_5 <- t.test(type_ddi2,type_ddi5)
t_2_5 |> glance()
t_2_5 |> report() #p < .001
#t_2_6 <- t.test(g2,g6)
#t_2_6 |> glance()
#t_2_6 |> report() # Not enough observations
t_3_4 <- t.test(type_ddi3,type_ddi4)
t_3_4 |> glance()
t_3_4 |> report() #p < .001
t_3_5 <- t.test(type_ddi3,type_ddi5)
t_3_5 |> glance()
t_3_5 |> report() #p < .001
#t_3_6 <- t.test(type_ddi3,type_ddi6)
#t_3_6 |> glance()
#t_3_6 |> report() # Not enough observations
t_4_5 <- t.test(type_ddi4,type_ddi5)
t_4_5 |> glance()
t_4_5 |> report() #p < .001
#t_4_6 <- t.test(type_ddi4,type_ddi6)
#t_4_6 |> glance()
#t_4_6 |> report() # Not enough observations
#t_5_6 <- t.test(type_ddi5,type_ddi6)
#t_5_6 |> glance()
#t_5_6 |> report() # Not enough observations

r1_2 <- c(type_ddi[1], type_ddi[2], "Welch Two Sample t-test", round(t_1_2$statistic, 3), round(t_1_2$p.value, 3), round(t_1_2$parameter, 3))
r1_3 <- c(type_ddi[1], type_ddi[3], "Welch Two Sample t-test", round(t_1_3$statistic, 3), round(t_1_3$p.value, 3), round(t_1_3$parameter, 3))
r1_4 <- c(type_ddi[1], type_ddi[4], "Welch Two Sample t-test", round(t_1_4$statistic, 3), round(t_1_4$p.value, 3), round(t_1_4$parameter, 3))
r1_5 <- c(type_ddi[1], type_ddi[5], "Welch Two Sample t-test", round(t_1_5$statistic, 3), round(t_1_5$p.value, 3), round(t_1_5$parameter, 3))
r1_6 <- c(type_ddi[1], type_ddi[6], "Not enough observations", "NA", "NA", "NA")
r2_3 <- c(type_ddi[2], type_ddi[3], "Welch Two Sample t-test", round(t_2_3$statistic, 3), round(t_2_3$p.value, 3), round(t_2_3$parameter, 3))
r2_4 <- c(type_ddi[2], type_ddi[4], "Welch Two Sample t-test", round(t_2_4$statistic, 3), round(t_2_4$p.value, 3), round(t_2_4$parameter, 3))
r2_5 <- c(type_ddi[2], type_ddi[5], "Welch Two Sample t-test", round(t_2_5$statistic, 3), round(t_2_5$p.value, 3), round(t_2_5$parameter, 3))
r2_6 <- c(type_ddi[2], type_ddi[6], "Not enough observations", "NA", "NA", "NA")
r3_4 <- c(type_ddi[3], type_ddi[4], "Welch Two Sample t-test", round(t_3_4$statistic, 3), round(t_3_4$p.value, 3), round(t_3_4$parameter, 3))
r3_5 <- c(type_ddi[3], type_ddi[5], "Welch Two Sample t-test", round(t_3_5$statistic, 3), round(t_3_5$p.value, 3), round(t_3_5$parameter, 3))
r3_6 <- c(type_ddi[3], type_ddi[6], "Not enough observations", "NA", "NA", "NA")
r4_5 <- c(type_ddi[4], type_ddi[5], "Welch Two Sample t-test", round(t_4_5$statistic, 3), round(t_4_5$p.value, 3), round(t_4_5$parameter, 3))
r4_6 <- c(type_ddi[4], type_ddi[6], "Not enough observations", "NA", "NA", "NA")
r5_6 <- c(type_ddi[5], type_ddi[6], "Not enough observations", "NA", "NA", "NA")

sup.tab.s2 <- rbind(r1_2,r1_3,r1_4,r1_5,r1_6,
      r2_3,r2_4,r2_5,r2_6,
      r3_4,r3_5,r3_6,
      r4_5,r4_6,
      r5_6)
colnames(sup.tab.s2) <- c("Group 1","Group 2","Statistical test","Statistic","p-value","df")
sup.tab.s2 |> write.csv(file = paste0("results/",set_name,"/sup.tab.s2.fig4a",set_name,".csv"),row.names = FALSE)

fig.4a.df1 <- fig.4a.df|> 
  group_by(type) |> 
  summarise(Mean.rate=mean(sigma.rate,na.rm=TRUE),
            sem.rate=std.error(sigma.rate,na.rm=TRUE),
            count=n()) |> 
  ungroup() |> 
  mutate(type=as.factor(type)) 

fig.4a <- fig.4a.df1 |> 
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
  ylim(0,0.0020)+
  scale_x_discrete(limits=rev)

antagonism_score <- -1
synergy_score <- +1

fig.4b.df <- df |> 
  select(sum_g,sigma.rate) |> 
  filter(!is.na(sum_g)) |> 
  group_by(sum_g) |> 
  summarise(Mean.rate=mean(sigma.rate,na.rm=TRUE),
            sem.rate=std.error(sigma.rate,na.rm=TRUE),
            count=n()) |> 
  ungroup() 

fig.4b <- fig.4b.df |> 
  ggplot(aes(x=sum_g,y=Mean.rate))+
  geom_bar(position=position_dodge(), stat="identity",
           fill='#6191F2',width=0.5)+
  geom_errorbar(position=position_dodge(0.5),aes(ymin=Mean.rate-sem.rate, ymax=Mean.rate+sem.rate), width=.2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+
  xlab("Score (More synergistic <==> More antagonistic)")+
  ylab(expression("Sigma rate ("~Bliss^2/MYA~")"))+
  theme(legend.position = "none")+
  geom_text(aes(label = paste0("n=",count)),vjust=-0.75)

# Group fig 4b
fig.4b.df <- df |> 
  select(sum_g,sigma.rate) |> 
  filter(!is.na(sum_g)) |> 
  mutate(sum_g1=cut(sum_g,breaks=c(-7,-3,+2,+7))) 

fig.4b.df1 <- fig.4b.df |> 
  group_by(sum_g1) |> 
  summarise(Mean.rate=mean(sigma.rate),
            sem.rate=std.error(sigma.rate),
            count=n()) |> 
  ungroup() 

s <- fig.4b.df |> filter(sum_g1=="(2,7]") |> select(sigma.rate) |> pull()
ad <- fig.4b.df |> filter(sum_g1=="(-3,2]") |> select(sigma.rate) |> pull()
an <- fig.4b.df |> filter(sum_g1=="(-7,-3]") |> select(sigma.rate) |> pull()

t1 <- t.test(s,ad) #p < .001
t1 |> glance()
t2 <- t.test(s,an) #p < .001
t2 |> glance()
t3 <- t.test(ad,an) #p = 0.001
t3 |> glance()

r1 <- c("Synergy", "Additivity", "Welch Two Sample t-test", round(t1$statistic, 3), round(t1$p.value, 3), round(t1$parameter, 3))
r2 <- c("Synergy", "Antagonism", "Welch Two Sample t-test", round(t2$statistic, 3), round(t2$p.value, 3), round(t2$parameter, 3))
r3 <- c("Additivity", "Antagonism", "Welch Two Sample t-test", round(t3$statistic, 3), round(t3$p.value, 3), round(t3$parameter, 3))


sup.tab.fig4b <- rbind(r1,r2,r3)
colnames(sup.tab.fig4b) <- c("Group 1","Group 2","Statistical test","Statistic","p-value","df")
sup.tab.fig4b |> write.csv(file = paste0("results/",set_name,"/sup.tab.s2.fig4b",set_name,".csv"),row.names = FALSE)

fig.4b.2 <- fig.4b.df1 |> 
  ggplot(aes(x=sum_g1,y=Mean.rate))+
  geom_bar(position=position_dodge(), stat="identity",
           fill='#6191F2',width=0.5)+
  geom_errorbar(position=position_dodge(0.5),aes(ymin=Mean.rate-sem.rate, ymax=Mean.rate+sem.rate), width=.2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+
  xlab("Score (More synergistic <==> More antagonistic)")+
  ylab(expression("Sigma rate ("~Bliss^2/MYA~")"))+
  theme(legend.position = "none")+
  geom_text(aes(label = paste0("n=",count)),vjust=-0.75)

# Rates ~ Interaction Type (ebw) (PPI). 5C
fig.4c.df <- df_target_DDI |>
  select(network,int_sign_ebw.DDI,starts_with("sigma.rate")) |> 
  select(network,int_sign_ebw.DDI,ends_with(paste0(set_col,".DDI"))) |> 
  rename(sigma.rate=starts_with("sigma.rate")) |> 
  filter(!is.na(int_sign_ebw.DDI)) 

fig.4c <- fig.4c.df |>
  ggline(x = "int_sign_ebw.DDI", y = "sigma.rate", add = "mean_se",
         group = "network", col="network")+
  xlab("Type of interaction (ebw)")+
  ylab(expression("Sigma rate ("~Bliss^2/MYA~")"))+
  theme_minimal()+
  theme(legend.position = "bottom")

fig.4c.df.n1 <- fig.4c.df |> 
  filter(network=="Co-functional (EcoCyc/GO-BP)") 

n1.s <- fig.4c.df.n1 |> filter(int_sign_ebw.DDI=="Synergy") |> select(sigma.rate) |> pull()
n1.ad <- fig.4c.df.n1 |> filter(int_sign_ebw.DDI=="Additivity") |> select(sigma.rate) |> pull()
n1.an <- fig.4c.df.n1 |> filter(int_sign_ebw.DDI=="Antagonism") |> select(sigma.rate) |> pull()

t1.n1 <- t.test(n1.s,n1.ad) #p < .001
t1.n1 |> glance()
t2.n1 <- t.test(n1.s,n1.an) #p < .001
t2.n1 |> glance()
t3.n1 <- t.test(n1.ad,n1.an) #p = 0.001775
t3.n1 |> glance()

fig.4c.df.n2 <- fig.4c.df |> 
  filter(network=="Small/medium-scale PPI")
n2.s <- fig.4c.df.n2 |> filter(int_sign_ebw.DDI=="Synergy") |> select(sigma.rate) |> pull()
n2.ad <- fig.4c.df.n2 |> filter(int_sign_ebw.DDI=="Additivity") |> select(sigma.rate) |> pull()
n2.an <- fig.4c.df.n2 |> filter(int_sign_ebw.DDI=="Antagonism") |> select(sigma.rate) |> pull()

t1.n2 <- t.test(n2.s,n2.ad) #p < .001
t1.n2 |> glance()
t2.n2 <- t.test(n2.s,n2.an) #p < .001
t2.n2 |> glance()
t3.n2 <- t.test(n2.ad,n2.an) #p = 0.852
t3.n2 |> glance()

r1 <- c("Synergy", "Additivity", "Welch Two Sample t-test", round(t1.n1$statistic, 3), round(t1.n1$p.value, 3), round(t1.n1$parameter, 3),"Co-functional (EcoCyc/GO-BP)")
r2 <- c("Synergy", "Antagonism", "Welch Two Sample t-test", round(t2.n1$statistic, 3), round(t2.n1$p.value, 3), round(t2.n1$parameter, 3),"Co-functional (EcoCyc/GO-BP)")
r3 <- c("Additivity", "Antagonism", "Welch Two Sample t-test", round(t3.n1$statistic, 3), round(t3.n1$p.value, 3), round(t3.n1$parameter, 3),"Co-functional (EcoCyc/GO-BP)")
r4 <- c("Synergy", "Additivity", "Welch Two Sample t-test", round(t1.n2$statistic, 3), round(t1.n2$p.value, 3), round(t1.n2$parameter, 3),"Small/medium-scale PPI")
r5 <- c("Synergy", "Antagonism", "Welch Two Sample t-test", round(t2.n2$statistic, 3), round(t2.n2$p.value, 3), round(t2.n2$parameter, 3),"Small/medium-scale PPI")
r6 <- c("Additivity", "Antagonism", "Welch Two Sample t-test", round(t3.n2$statistic, 3), round(t3.n2$p.value, 3), round(t3.n2$parameter, 3),"Small/medium-scale PPI")

sup.tab.fig4c <- rbind(r1,r2,r3,r4,r5,r6)
colnames(sup.tab.fig4c) <- c("Group 1","Group 2","Statistical test","Statistic","p-value","df","Network")
sup.tab.fig4c |> write.csv(file = paste0("results/",set_name,"/sup.tab.s2.fig4c",set_name,".csv"),row.names = FALSE)

#Red:metabolic,blue:PPI
# Rates ~ path length. 5D
fig.4d.df <- df_target_DDI |>
  select(network,path.length,starts_with("sigma.rate")) |> 
  select(network,path.length,ends_with(paste0(set_col,".DDI"))) |> 
  rename(sigma.rate=starts_with("sigma.rate")) |> 
  filter(!is.na(path.length)) |> 
  filter(network==net[1]|network==net[3]) |> 
  mutate(path.length.corr=ifelse(is.infinite(path.length),"Not connected",
                                 ifelse(path.length==0,"0",
                                        ifelse(path.length==1,"1",
                                               ifelse(path.length==2,"2",
                                                      ifelse(path.length==3,"3",
                                                             ifelse(path.length==4,"4",
                                                                    ifelse(path.length==5,"5",
                                                                           ifelse(path.length>=6,">6","Error"))))))))) |> 
  mutate(path.length.corr=factor(path.length.corr,levels=c("0","1","2","3","4","5",">6","Not connected","Error"))) 

fig.4d.df1 <- fig.4d.df |> 
  filter(path.length.corr!="Not connected") |> 
  mutate(path.length.corr1=ifelse(path.length.corr==">6",6,ifelse(path.length.corr=="5",5,ifelse(path.length.corr=="4",4,ifelse(path.length.corr=="3",3,ifelse(path.length.corr=="2",2,ifelse(path.length.corr=="1",1,"error"))))))) |> 
  select(network,sigma.rate,path.length.corr1) |> mutate(path.length.corr1=as.numeric(path.length.corr1))

fig.4d.df1.n1 <- fig.4d.df1 |> filter(network=="Co-functional (EcoCyc/GO-BP)") |> group_by(path.length.corr1) |> summarise(Mean.sigma=mean(sigma.rate,na.rm=T),
                                                                                                                           sem.sigma=std.error(sigma.rate),
                                                                                                                           count=paste0("n=",n()))
mod1 <- lm(Mean.sigma~path.length.corr1,fig.4d.df1.n1)
mod1 |> tidy()
mod1 |> glance()
mod1 |> report()
mod1.t <- mod1 |> tidy()

fig.4d.df1.n2 <- fig.4d.df1 |> filter(network=="Small/medium-scale PPI") |> group_by(path.length.corr1) |> summarise(Mean.sigma=mean(sigma.rate,na.rm=T),
                                                                                                                     sem.sigma=std.error(sigma.rate),
                                                                                                                     count=paste0("n=",n()))
mod2 <- lm(Mean.sigma~path.length.corr1,fig.4d.df1.n2)
mod2 |> tidy()
mod2 |> glance()
mod2 |> report()
mod2.t <- mod2 |> tidy()


fig.4d <- fig.4d.df |>
  ggline(x = "path.length.corr", y = "sigma.rate", add = c("mean_se"),
         group = "network", col="network",linetype = 1)+
  xlab("Minimum path length between connected targets")+
  ylab(expression("Sigma rate ("~Bliss^2/MYA~")"))+
  theme_minimal()+
  theme(legend.position = "none")+ 
  scale_x_discrete(labels = c("1","2","3","4","5", "\u2265 6","Not connected"))

r1 <- c("OLS(Mean.sigma ~ path.length.corr", "Intercept", as.character(c(round(mod1.t[1,2], 5), round(mod1.t[1,3], 6), round(mod1.t[1,4], 2),round(mod1.t[1,5], 5))),"Co-functional (EcoCyc/GO-BP)")
r2 <- c("OLS(Mean.sigma ~ path.length.corr", "Slope", as.character(c(round(mod1.t[2,2], 5), round(mod1.t[2,3], 6), round(mod1.t[2,4], 2),round(mod1.t[2,5], 5))),"Co-functional (EcoCyc/GO-BP)")
r3 <- c("OLS(Mean.sigma ~ path.length.corr", "Intercept", as.character(c(round(mod2.t[1,2], 5), round(mod2.t[1,3], 6), round(mod2.t[1,4], 2),round(mod2.t[1,5], 5))),"Small/medium-scale PPI")
r4 <- c("OLS(Mean.sigma ~ path.length.corr", "Slope", as.character(c(round(mod2.t[2,2], 5), round(mod2.t[2,3], 6), round(mod2.t[2,4], 2),round(mod2.t[2,5], 5))),"Small/medium-scale PPI")

sup.tab.fig4d <- rbind(r1,r2,r3,r4)
colnames(sup.tab.fig4d) <- c("Model","Term","Estimate","std.error","statistic","value","Network")
sup.tab.fig4d |> write.csv(file = paste0("results/",set_name,"/sup.tab.s2.fig4d",set_name,".csv"),row.names = FALSE)

# Fig 4
fig4 <- wrap_plots(fig.4a,fig.4b,
                   fig.4c,fig.4d,nrow=2)+ 
  plot_annotation(tag_levels = 'A',tag_suffix = '.')&
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 15, hjust = 0, vjust = 0))+
  theme(axis.title = element_text(size = 12),  # Set axis title size
        axis.text = element_text(size = 11))

ggsave(
  filename = paste0("results/", set_name, "/fig4_", set_name, ".pdf"),
  plot = fig4,
  width = 15,
  height = 15,
  dpi=600
)

ggsave(
  filename = paste0("results/", set_name, "/fig4_", set_name, ".png"),
  plot = fig4,
  width = 15,
  height = 15,
  dpi=600
)

# Rates ~ Interaction Type (all strains) (all networks).
df_target_DDI |> 
  select(network,sum_g.target,starts_with("sigma.rate")) |> 
  select(network,sum_g.target,ends_with(paste0(set_col,".DDI"))) |> 
  rename(sigma.rate=starts_with("sigma.rate")) |> 
  filter(!is.na(sigma.rate)) |> 
  ggline(x = "sum_g.target", y = "sigma.rate", add = "mean_se",
         group = "network", palette = "jco")+
  stat_compare_means(aes(group = network), label = "p.signif",
                     label.y = c(40, 40, 40))+
  facet_wrap(~network)+
  xlab("Type of interaction")+
  ylab("Sigma rate")+
  theme_minimal()

# Mean.path. Sup_D
sup_d <- df_DDI_tot |>
  filter(network==net[1]|network==net[3]) |> 
  filter(!is.infinite(mean.path.length)) |> 
  mutate(type=ifelse(is.na(type),"Non-targets",type)) |>
  select(network,DRUG_ID,sum_g,connection_onoff,type,mean.path.length,
         starts_with(c("sigma.rate","clusters"))) |> 
  select(network,DRUG_ID,connection_onoff,type,mean.path.length,
         ends_with(paste0(set_col,".DDI"))) |> 
  rename(sigma.rate=starts_with("sigma.rate"),
         clusters=starts_with("clusters")) |> 
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
  ylab("Minimum path distance between targets")+
  theme(legend.position = "none")+
  geom_text(aes(label = count,vjust=-0.75,hjust=1))+
  ylim(0,6)

#Mean.K.edge
sup_e <- df_DDI_tot |>
  filter(network==net[1]|network==net[3]) |> 
  mutate(type=ifelse(is.na(type),"Non-targets",type)) |>
  select(network,DRUG_ID,connection_onoff,type,mean.path.length,mean.k.edge,
         starts_with(c("sigma.rate","clusters"))) |> 
  select(network,DRUG_ID,connection_onoff,type,mean.path.length,mean.k.edge,
         ends_with(paste0(set_col,".DDI"))) |> 
  rename(sigma.rate=starts_with("sigma.rate"),
         clusters=starts_with("clusters")) |>
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
  ylab("Average K-edge connectivity between targets")+
  theme(legend.position = "bottom")+
  geom_text(aes(label = count,vjust=-0.75,hjust=1))+
  ylim(0,35)

# Mean Deg
sup_f <- df_DDI_tot |>
  filter(network==net[1]|network==net[3]) |> 
  mutate(type=ifelse(is.na(type),"Non-targets",type)) |>
  select(network,DRUG_ID,connection_onoff,type,mean.path.length,mean.mean.degree,
         starts_with(c("sigma.rate","clusters"))) |> 
  select(network,DRUG_ID,connection_onoff,type,mean.path.length,mean.mean.degree,
         ends_with(paste0(set_col,".DDI"))) |> 
  rename(sigma.rate=starts_with("sigma.rate"),
         clusters=starts_with("clusters")) |> 
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
  ylab("Average node degree")+
  theme(legend.position = "none")+
  geom_text(aes(label = count,vjust=-0.75,hjust=1))+
  ylim(0,45)

# Mean Betweenness centrality
sup_g <- df_DDI_tot |>
  filter(network==net[1]|network==net[3]) |> 
  mutate(type=ifelse(is.na(type),"Non-targets",type)) |>
  select(network,DRUG_ID,connection_onoff,type,mean.mean.bet,mean.mean.degree,
         starts_with(c("sigma.rate","clusters"))) |> 
  select(network,DRUG_ID,connection_onoff,type,mean.mean.bet,mean.mean.degree,
         ends_with(paste0(set_col,".DDI"))) |> 
  rename(sigma.rate=starts_with("sigma.rate"),
         clusters=starts_with("clusters")) |> 
  distinct() |> 
  group_by(type,network) |> 
  summarise(Mean.mean_bet=mean(mean.mean.bet,na.rm=T),
            sem.mean_bet=std.error(mean.mean.bet),
            count=paste0("n=",n())) |> 
  ungroup() |> 
  ggplot(aes(x=type,y=Mean.mean_bet,fill=network))+
  geom_bar(position=position_dodge(), stat="identity",
           colour='black',width=0.5)+
  geom_errorbar(position=position_dodge(0.5),aes(ymin=Mean.mean_bet-sem.mean_bet, ymax=Mean.mean_bet+sem.mean_bet), width=.2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+
  xlab("Interaction type across all strains")+
  ylab("Average betweenness centrality")+
  theme(legend.position = "none")+
  geom_text(aes(label = count,vjust=-0.75,hjust=1))

# Mean Eigenvector centrality
sup_h <- df_DDI_tot |>
  filter(network==net[1]|network==net[3]) |> 
  mutate(type=ifelse(is.na(type),"Non-targets",type)) |>
  select(network,DRUG_ID,connection_onoff,type,mean.path.length,mean.mean.EVC,
         starts_with(c("sigma.rate","clusters"))) |> 
  select(network,DRUG_ID,connection_onoff,type,mean.path.length,mean.mean.EVC,
         ends_with(paste0(set_col,".DDI"))) |> 
  rename(sigma.rate=starts_with("sigma.rate"),
         clusters=starts_with("clusters")) |> 
  distinct() |> 
  group_by(type,network) |> 
  summarise(Mean.mean_EVC=mean(mean.mean.EVC,na.rm=T),
            sem.mean_EVC=std.error(mean.mean.EVC),
            count=paste0("n=",n())) |> 
  ungroup() |> 
  ggplot(aes(x=type,y=Mean.mean_EVC,fill=network))+
  geom_bar(position=position_dodge(), stat="identity",
           colour='black',width=0.5)+
  geom_errorbar(position=position_dodge(0.5),aes(ymin=Mean.mean_EVC-sem.mean_EVC, ymax=Mean.mean_EVC+sem.mean_EVC), width=.2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+
  xlab("Interaction type across all strains")+
  ylab("Average eigenvector centrality")+
  theme(legend.position = "none")+
  geom_text(aes(label = count,vjust=-0.75,hjust=1))

sup_k <- df_DDI_tot |>
  filter(network==net[1]|network==net[3]) |> 
  mutate(type=ifelse(is.na(type),"Non-targets",type)) |>
  select(network,DRUG_ID,connection_onoff,type,mean.path.length,mean.mean.degree,mean.k.edge,starts_with(c("sigma.rate","clusters"))) |> 
  select(network,DRUG_ID,connection_onoff,type,mean.path.length,mean.mean.degree,mean.k.edge, ends_with(paste0(set_col))) |> 
  rename(sigma.rate=starts_with("sigma.rate"),
         clusters=starts_with("clusters"))|> 
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
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))

sup_l <- df_DDI_tot |>
  filter(network==net[1]|network==net[3]) |> 
  mutate(type=ifelse(is.na(type),"Non-targets",type)) |>
  select(network,DRUG_ID,connection_onoff,type,mean.path.length,mean.mean.degree,mean.k.edge,starts_with(c("sigma.rate","clusters"))) |> 
  select(network,DRUG_ID,connection_onoff,type,mean.path.length,mean.mean.degree,mean.k.edge, ends_with(paste0(set_col))) |> 
  rename(sigma.rate=starts_with("sigma.rate"),
         clusters=starts_with("clusters")) |> 
  distinct() |> 
  mutate(mean.deg.group=ifelse(mean.mean.degree<=5,"0-5",
                               ifelse(mean.mean.degree<=10,"5-10",
                                      ifelse(mean.mean.degree<=20,"10-20",
                                             ifelse(mean.mean.degree>=20,">20","error"))))) |> 
  mutate(mean.deg.group=factor(mean.deg.group,levels=c("0-5","5-10","10-20",">20","error"))) |>
  ggline(x = "mean.deg.group", y = "sigma.rate", add = "mean_se",
         group = "network", col="network")+
  xlab("Average node degree")+
  ylab(expression("Sigma rate ("~Bliss^2/MYA~")"))+
  theme_minimal()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))

sup_m <- df_DDI_tot |>
  filter(network==net[1]|network==net[3]) |> 
  mutate(type=ifelse(is.na(type),"Non-targets",type)) |>
  select(network,DRUG_ID,connection_onoff,type,mean.path.length,mean.mean.bet,mean.k.edge,starts_with(c("sigma.rate","clusters"))) |> 
  select(network,DRUG_ID,connection_onoff,type,mean.path.length,mean.mean.bet,mean.k.edge, ends_with(paste0(set_col))) |> 
  rename(sigma.rate=starts_with("sigma.rate"),
         clusters=starts_with("clusters")) |> 
  distinct() |> 
  mutate(mean.bet.group=ifelse(mean.mean.bet<=100,"0-100",
                               ifelse(mean.mean.bet<=250,"100-250",
                                      ifelse(mean.mean.bet<=500,"250-500",
                                             ifelse(mean.mean.bet<=1000,"500-1000",
                                                    ifelse(mean.mean.bet<=2000,"1000-2000",
                                                           ifelse(mean.mean.bet<=3000,"2000-3000",
                                                                  ifelse(mean.mean.bet>=3000,">3000","error")))))))) |> 
  mutate(mean.bet.group=factor(mean.bet.group,levels=c("0-100","100-250",
                                                       "250-500","500-1000","1000-2000","2000-3000",">3000","error"))) |>
  ggline(x = "mean.bet.group", y = "sigma.rate", add = "mean_se",
         group = "network", col="network")+
  xlab("Average betweenness centrality")+
  ylab(expression("Sigma rate ("~Bliss^2/MYA~")"))+
  theme_minimal()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))

sup_n <- df_DDI_tot |>
  filter(network==net[1]|network==net[3]) |> 
  mutate(type=ifelse(is.na(type),"Non-targets",type)) |>
  select(network,DRUG_ID,connection_onoff,type,mean.path.length,mean.mean.EVC,mean.k.edge,starts_with(c("sigma.rate","clusters"))) |> 
  select(network,DRUG_ID,connection_onoff,type,mean.path.length,mean.mean.EVC,mean.k.edge, ends_with(paste0(set_col))) |> 
  rename(sigma.rate=starts_with("sigma.rate"),
         clusters=starts_with("clusters")) |> 
  distinct() |> 
  mutate(mean.EVC.group=ifelse(mean.mean.EVC<=0.025,"0-0.025",
                               ifelse(mean.mean.EVC<=0.05,"0.025-0.05",
                                      ifelse(mean.mean.EVC<=0.1,"0.05-0.1",
                                             ifelse(mean.mean.EVC<=0.2,"0.1-0.2",
                                                    ifelse(mean.mean.EVC<=0.3,"0.2-0.3",
                                                           ifelse(mean.mean.EVC>=0.3,">0.3","error"))))))) |> 
  mutate(mean.EVC.group=factor(mean.EVC.group,levels=c("0-0.025","0.025-0.05","0.05-0.1","0.1-0.2","0.2-0.3",">0.3","error"))) |>
  ggline(x = "mean.EVC.group", y = "sigma.rate", add = "mean_se",
         group = "network", col="network")+
  xlab("Average eigenvector centrality")+
  ylab(expression("Sigma rate ("~Bliss^2/MYA~")"))+
  theme_minimal()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))

# Make figures:
# 20 x 20. SupFig
SupFig <- wrap_plots(sup_a,sup_b,sup_c,
                     sup_d,sup_e,sup_f,
                     sup_g,sup_h,sup_i,
                     sup_j,sup_k,sup_l,
                     sup_m,sup_n,
                     ncol=3)+ 
  plot_annotation(tag_levels = 'A',tag_suffix = '.')&
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 15, hjust = 0, vjust = 0)) 

df_DDI_tot |> 
  filter(network==net[1]|network==net[3]) |> 
  filter(Drug1 %in% c("A22","Novobiocin","Ciprofloxacin",
                      "Meropenem","Sulfamonomethoxine","Trimethoprim")) |>
  filter(Drug2 %in% c("A22","Novobiocin","Ciprofloxacin",
                      "Meropenem","Sulfamonomethoxine","Trimethoprim")) |> 
  select(network,DRUG_ID,clusters_all,sigma.rate_all,clusters_set,sigma.rate_set,
         mean.path.length,mean.k.edge,mean.mean.degree,
         starts_with("int"),starts_with("SR")) |> 
  filter(DRUG_ID %in% c("Sulfamonomethoxine-Trimethoprim",
                        "A22-Novobiocin",
                        "Ciprofloxacin-Meropenem")) |> 
  arrange(DRUG_ID,network) |> 
  t() |> 
  write.csv(paste0("results/",set_name,"/sup_table1_",set_col,".csv"),col.names = FALSE)

#Rates ~ interaction type and resistance type.
df_DDI_tot |> ggline(x = "int_sign_ebw", y = "sigma.rate_all", add = "mean_se",
                     group = "network", palette = "jco")+
  stat_compare_means(aes(group = network), label = "p.signif",
                     label.y = c(40, 40, 40))+
  facet_wrap(~network)+
  xlab("Type of interaction")+
  ylab("Sigma rate")+
  theme_minimal()


df1 <- df_DDI_tot |> 
  filter(network==net[1]|network==net[3]) |> 
  filter(DRUG_ID %in% c("A22-Novobiocin","Ciprofloxacin-Meropenem",
                        "Sulfamonomethoxine-Trimethoprim","Doxycycline-Trimethoprim"))

df1 |> t() |> write.csv("results/allddi/example_DDI.csv")




a <- 8
ggsave(
  filename = paste0("results/", set_name, "/supfig_1", set_name, ".pdf"),
  plot = SupFig,
  width = 2 * a,
  height = 3 * a,
  dpi=600
)

ggsave(
  filename = paste0("results/", set_name, "/supfig_1", set_name, ".png"),
  plot = SupFig,
  width = 2 * a,
  height = 3 * a,
  dpi=600
)

