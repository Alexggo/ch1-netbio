#!/usr/bin/env Rscript

## ----load---------------------------------------------------------------------
library(pacman)
p_load(igraph,tidyverse,ape,geomorph,matrixStats,
       ggpubr,plotrix,tidymodels,patchwork,here)


## -----------------------------------------------------------------------------
filenames <- c("EcoCyc.goldstandardset.txt", "EcoliNet.v1.txt",
                       "GN.INT.EcoliNet.3568gene.23439link.txt",
                       "GO-BP.goldstandardset.txt",
                       "CC.EcoliNet.v1.2296gene.50528link.txt",
                       "CX.INT.EcoliNet.v1.4039gene.67494link.txt",
                       "DC.EcoliNet.2283gene.9643link.txt",
                       "EcoliNet.v1.benchmark.txt",
                       "HT.INT.EcoliNet.3209gene.15543link.txt",
                       "LC.INT.EcoliNet.764gene.1073link.txt",
                       "PG.INT.EcoliNet.v1.1817gene.17504link.txt")

# Ecocyc_goldstandard:1, Small-medium PPI:10
filepath_ecolinet <- file.path("data/5.Targets_NetworkDistance", filenames)
net <- c("Co-functional (EcoCyc)", "EN: all networks",
         "Similar genomic context", "Co-functional (GO-BP)",
         "Co-citation", "Co-expression",
         "Co-occurence of prot. domains", "Co-functional (EcoCyc/GO-BP)",
         "High-throughput PPI", "Small/medium-scale PPI",
         "Similar phylogenetic profiles")

list_ecolinet <- list()
for (i in seq_along(filepath_ecolinet)) {
  list_ecolinet[[i]] <- read_table(filepath_ecolinet[[i]])
  list_ecolinet[[i]]$node1 <- paste0("eco:", list_ecolinet[[i]]$node1)
  list_ecolinet[[i]]$node2 <- paste0("eco:", list_ecolinet[[i]]$node2)
}

list_net <- list_ecolinet

list_net <- list_net |>
  map(select,node1,node2) 

names(list_net) <- net


## -----------------------------------------------------------------------------
list_graphs <- list_net |> 
  map(graph_from_data_frame, directed=FALSE, vertices=NULL)

l_deg_list <- list_graphs |> 
  map(igraph::degree) 

l_deg <- list_graphs |> 
  map_dfr(igraph::degree) |> 
  t() |> 
  as.data.frame() 

l_deg <- tibble::rownames_to_column(l_deg, "Node")
colnames(l_deg) <- c("Node", filenames)

l_deg_long <- l_deg |> 
  pivot_longer(names_to = "Network", values_to = "Degree", 2:12) |> 
  drop_na()

before_trimming <- l_deg_long |>
  ggplot(aes(x=Degree,y=Network,fill=Network)) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_y_discrete(limits=rev)

# Trim outliers
# q_table <- l_deg_long |>
#   group_by(Network) |>
#   summarise(P25=quantile(Degree, probs = 0.25),
#             P50=quantile(Degree, probs = 0.50),
#             P75=quantile(Degree, probs = 0.75),
#             P80=quantile(Degree, probs = 0.8),
#             P85=quantile(Degree, probs = 0.85),
#             P90=quantile(Degree, probs = 0.9),
#             P95=quantile(Degree, probs = 0.95),
#             P975=quantile(Degree, probs = 0.975),
#             P100=quantile(Degree, probs = 1))
# 
# 
# a <- l_deg_list |> 
#   map(as.data.frame) |> 
#   map(rownames_to_column, "Node") |> 
#   map(rename, "Degree"=`.x[[i]]`)
# for (i in seq_along(filenames)){
#   a[[i]] <- a[[i]] |> filter(Degree<q_table$P975[i]) 
# }
# names(a) <- filenames
# l_deg_long_trimmed <- bind_rows(a,.id="Network") |> 
#   drop_na()
# 
# after_trimming <- l_deg_long_trimmed |> 
#   ggplot(aes(x=Degree,y=Network,fill=Network))+
#   geom_boxplot()+
#   theme_minimal()+
#   theme(legend.position = "none")+
#   scale_y_discrete(limits=rev)
# 
# wrap_plots(before_trimming,after_trimming)
# 
# #Remove nodes with degree greater than 90%.
# list_graphs <- list_net |> 
#   map(graph_from_data_frame,directed=FALSE,vertices=NULL)
# trimmed_graph <- list()
# for (i in 1:length(l_deg_list)){
#   vect <- l_deg_list[[i]] |> sort()
#   from <- round(1*length(vect))
#   to <- length(vect)
#   print(paste("From",from,"to",to))
#   trimmed <- l_deg_list[[i]][from:to] |> 
#     as.data.frame()
#   trimmed <- rownames_to_column(trimmed, "Node")
#   colnames(trimmed) <- c("Node","Degree")
#   trimmed <- trimmed |> arrange(Degree)
#   trimmed_graph[[i]] <- delete_vertices(list_graphs[[i]],trimmed$Node)
# }
# list_graphs <- trimmed_graph

l_deg_list <- list_graphs |> 
  map(igraph::degree) 

# Which proteins are targeted by drugs?
nodes <- read.csv("data/5.Targets_NetworkDistance/DrugTargets3_ecoli.csv") |> 
  select(Drug,KEGG_eco) |>
  distinct() |>
  arrange(Drug) |> 
  filter(KEGG_eco!="")
possible_targets <- nodes$KEGG_eco |> unique()

### PLOTTING NETWORKS.
#Before trimming:
# Nodes:885,4099,3568,1474,2296,4039,2283,1835,3209,764,1817,1165,1157
# Targets:10,25,25,17,12,25,6,20,25,18,15,9,9
# Drugs:16,36,36,31,23,36,16,31,36,31,26,17,17

# After P90% trimming:
# Nodes: 795,3688,3210,1326,2065,3634,2054,1651,2887,687,1634,1047,1040
# Targets: 10,24,25,15,12,24,6,18,25,17,14,9,4
# Drugs: 16,34,36,29,23,34,16,29,36,30,25,17,8

#After P85% trimming:
# Nodes: 751,3484,3032,1252,1951,3432,1940,1559,2727,648,1543,989,982
# Targets:10,24,25,15,12,23,5,18,24,17,12,9,4
# Drugs:16,34,36,29,23,34,14,29,35,30,24,17,8

# Export network for cytoscape
# 
# library(KeyPathwayMineR)
# igraph_to_sif(list_graphs[[1]], "data/5.Targets_NetworkDistance/ecocycnet.sif")


## -----------------------------------------------------------------------------
pdf(file = file.path("networks","network_graphs.pdf"))

for (i in 1:length(list_graphs)){
  g1 <- list_graphs[[i]]
  names <- V(g1) |> names()
  print(net[i])
  value1 <- names %in% possible_targets |> sum()
  value2 <- nodes$Drug[possible_targets%in%names] |> unique() |> length()
  V(g1)$color <- ifelse(names %in% possible_targets, "red", "lightblue")
  V(g1)$size <- ifelse(V(g1) %in% possible_targets, 40, 5)
  V(g1)$label.cex <- 0.3
  v.number <- V(g1) |> length()
  mean.length <- average.path.length(g1)
  plot(g1,vertex.label="",vertex.size=2)
  mytitle = paste(net[i])
  mysubtitle1 = paste0("Average path length =",round(mean.length,4),". Number of vertexes =",v.number)
  mysubtitle2 = paste0("Number of targets =",value1,". Number of drugs =",value2)
  mtext(side=3, line=3, at=-0.07, adj=0, cex=1, mytitle)
  mtext(side=3, line=2, at=-0.07, adj=0, cex=0.7, mysubtitle1)
  mtext(side=3, line=1, at=-0.07, adj=0, cex=0.7, mysubtitle2)
  print(paste("Nodes",v.number,"Targets",value1,"Drugs",value2))
}
dev.off()

pdf(file = file.path("networks","network_graphs_notext.pdf"))

for (i in 1:length(list_graphs)){
  g1 <- list_graphs[[i]]
  names <- V(g1) |> names()
  print(net[i])
  value1 <- names %in% possible_targets |> sum()
  value2 <- nodes$Drug[possible_targets%in%names] |> unique() |> length()
  V(g1)$color <- ifelse(names %in% possible_targets, "red", "lightblue")
  V(g1)$size <- ifelse(V(g1) %in% possible_targets, 40, 5)
  V(g1)$label.cex <- 0.3
  v.number <- V(g1) |> length()
  mean.length <- average.path.length(g1)
  plot(g1,vertex.label="",vertex.size=2)
}
dev.off()


## -----------------------------------------------------------------------------
# Calculate path.length, k-edge and node degree for these targets.
mean_pathlength <- list_graphs |> 
  map_dbl(average.path.length)
shortest_paths <-  list_graphs |> 
  map(shortest.paths)

list_path <- list()
target_list <- list()
drug_list <- list()
for (i in 1:length(shortest_paths)){
  dist_mat <- shortest_paths[[i]]
  # Filter columns and rows for known targets
  r <- rownames(dist_mat) %in% possible_targets
  c <- colnames(dist_mat) %in% possible_targets
  # How many targets have drugs in Brochado?
  targets <- rownames(dist_mat)[rownames(dist_mat) %in% possible_targets]
  target_list[[i]] <- targets
  drugs <- nodes |> 
    filter(KEGG_eco%in%targets) |>
    select(Drug) |> 
    unique() |> 
    pull() |> 
    sort()
  print(paste("Network:",filenames[[i]]))
  Tar <-   rownames(dist_mat) %in% possible_targets |> sum() 
  print(paste("Number of targets:",Tar))
  print(targets)
  D <-   drugs |> length()
  print(paste("Number of drugs:",D))
  drugs |> print()
  # For small matrix of targets from Brochado, use dist_mat[r,c]
  distnew2 <- dist_mat |> 
    as.data.frame()
  # Convert matrix to table
  distnew2 <- distnew2[ , order(names(distnew2))]
  distnew2$KEGG1 <- rownames(distnew2)
  distnew2 <- distnew2 |>
    arrange(KEGG1) |> 
    select(-KEGG1)
  distnew2[lower.tri(distnew2)] <- NA
  distnew2$KEGG1 <- colnames(distnew2)
  distnew2 <- distnew2 |>
    select(dim(distnew2)[2],1:(dim(distnew2)[2]-1))
  list_path[[i]] <- distnew2  |> 
    pivot_longer(names_to="KEGG2",values_to="path.length",2:(dim(distnew2)[2])) |> 
    filter(!is.na(path.length))
}


#k.edge.connectivity between known nodes that are in g1
list_Kedge <- list()
for (i in 1:length(list_graphs)){
  # For small set use:
  #selection <- possible_targets[possible_targets %in% rownames(shortest_paths[[i]])]
  
  selection <- rownames(shortest_paths[[i]])
  t <- matrix(nrow = length(selection),ncol=length(selection))
  for (j in 1:length(selection)){
    for (k in 1:length(selection)){
      if(j==k){t[j,k] <- NA}else{
        t[j,k] <- edge_connectivity(list_graphs[[i]],source = selection[j],target = selection[k])
        colnames(t) <- selection
        rownames(t) <- selection
        
        connnew <- t |> as.data.frame()
        connnew <- connnew[ , order(names(connnew))]
        connnew$KEGG1 <- rownames(connnew)
        connnew <- connnew |>
          arrange(KEGG1) |> 
          select(-KEGG1)
        connnew[lower.tri(connnew)] <- NA
        connnew$KEGG1 <- colnames(connnew)
        connnew1 <- connnew |>
          select(dim(connnew)[2],1:(dim(connnew)[2]-1))
        connnew2 <- connnew1  |> 
          pivot_longer(names_to="KEGG2",values_to="K.edge",2:(dim(connnew)[2]))|>
          filter(!is.na(K.edge))
        list_Kedge[[i]] <- connnew2
        print(paste("i",i,"j",j,"k",k))
      }
    }
  }
}

# Merge both path length and K edge connectivity.
list_path <- list_path |> 
  map(mutate,KEGG1_KEGG2=paste0(KEGG1,"_",KEGG2))
list_Kedge <- list_Kedge |> 
  map(mutate,KEGG1_KEGG2=paste0(KEGG1,"_",KEGG2))

dist_conn <- map2(list_path, list_Kedge, full_join, by = "KEGG1_KEGG2") |> 
  map(select,KEGG1.x,KEGG2.x,path.length,K.edge) |> 
  map(rename,KEGG1=KEGG1.x,
      KEGG2=KEGG2.x,
      path.length=path.length) |> 
  map(arrange,K.edge) |> 
  map(distinct) |> 
  map(filter,KEGG1!=KEGG2)


l_deg_df <- l_deg_list |> map(as.data.frame)

for (i in 1:length(l_deg_df)){
  KEGG <- l_deg_df[[i]] |> rownames()
  degree <- l_deg_df[[i]][,1]
  l_deg_df[[i]] <- data.frame("KEGG1"=KEGG,
                              "Degree1"=degree) |> 
    distinct()
}

tabnew2 <- map2(dist_conn, l_deg_df, left_join, by = "KEGG1") 
names(tabnew2) <- net
l_deg_df2 <- l_deg_df
l_deg_df2 <- l_deg_df2 |> map(rename,"KEGG2"=KEGG1,"Degree2"=Degree1)
dist_conn_deg <- map2(tabnew2, l_deg_df2, left_join, by = "KEGG2") |> 
  map(mutate,KEGG1_KEGG2=paste0(KEGG1,"_",KEGG2))

sum_deg <- dist_conn_deg |> 
  map(pivot_longer,names_to="category",values_to="degree",5:6) |> 
  map(group_by,KEGG1,KEGG2) |> 
  map(summarise,mean_deg=mean(degree,na.rm=TRUE),#maybe round means.
      min_deg=min(degree,na.rm=TRUE),
      max_deg=max(degree,na.rm=TRUE))  |> 
  map(mutate,KEGG1_KEGG2=paste0(KEGG1,"_",KEGG2)) |> 
  map(ungroup) |> 
  map(distinct)

tabnew2 <- map2(dist_conn_deg, sum_deg, left_join, by = "KEGG1_KEGG2") |> 
  map(select,-c(KEGG1.y,KEGG2.y)) |> 
  map(rename,KEGG1=KEGG1.x,KEGG2=KEGG2.x)

nodes1 <- nodes |>
  rename(Drug1=Drug,
         KEGG1=KEGG_eco)

df_join <- tabnew2 |> 
  map(full_join,nodes1,by="KEGG1")

nodes2 <- nodes |>
  rename(Drug2=Drug,
         KEGG2=KEGG_eco)
df_join1 <- df_join |> 
  map(full_join,nodes2,by="KEGG2") |> 
  map(mutate,KEGG1_KEGG2=paste0(KEGG1,"_",KEGG2)) |> 
  map(distinct) |> 
  map(drop_na)


df_join2 <- list()
for (j in 1:length(df_join1)){
  z <- df_join1[[j]]
  for (i in 1:dim(z)[1]){
  a <- z[i,]$Drug1
  b <- z[i,]$Drug2
  ab <- sort(c(a,b))
  z[i,]$Drug1 <- ab[1]
  z[i,]$Drug2 <- ab[2]
}
  df_join2[[j]] <- z
}

df_join2 <- df_join2 |> 
  map(mutate,drug_pair=paste0(Drug1,"_",Drug2)) |> 
  map(select,KEGG1_KEGG2,KEGG1,KEGG2,drug_pair,Drug1,Drug2,path.length,K.edge,Degree1,Degree2,
      mean_deg,min_deg,max_deg) |> 
  map(distinct) |> 
  map(filter,Drug1!=Drug2) |> 
  map(filter,Drug1<Drug2)


adjacency <- list_graphs |> 
  map(as_adjacency_matrix,type = "both",names = TRUE,
      sparse = igraph_opt("sparsematrices")) |> 
  map(as.matrix) |> 
  map(as.data.frame)

df_adja <- list()
for (i in 1:length(adjacency)){
  adj <- adjacency[[i]]
  distnew2 <- adj[,order(colnames(adj))]
  distnew2$KEGG1 <- rownames(distnew2)
  distnew2 <- distnew2 |>
    arrange(KEGG1) |> 
    select(-KEGG1)
  distnew2[lower.tri(distnew2)] <- NA
  distnew2$KEGG1 <- colnames(distnew2)
  distnew2 <- distnew2 |>
    select(dim(distnew2)[2],1:(dim(distnew2)[2]-1))
  distnew3 <- distnew2  |> 
    pivot_longer(names_to="KEGG2",values_to="adjacency",2:(dim(distnew2)[2]))|>
    filter(!is.na(adjacency)) |> 
    mutate(KEGG1_KEGG2=paste0(KEGG1,"_",KEGG2))
  df_adja[[i]] <- distnew3
}

df_adja <- df_adja |> 
  map(distinct) |> 
  map(filter,KEGG1<KEGG2)

dist_conn_deg_adj <- map2(df_join2,df_adja,left_join,by="KEGG1_KEGG2") |> 
  map(select,-c(KEGG1.y,KEGG2.y)) |> 
  map(rename,KEGG1=KEGG1.x,KEGG2=KEGG2.x)
names(dist_conn_deg_adj) <- net



## -----------------------------------------------------------------------------
DDI_dist_conn_deg_adj <- dist_conn_deg_adj |> 
  map(dplyr::ungroup) |> 
  map(dplyr::group_by,drug_pair,Drug1,Drug2) |> 
  map(dplyr::summarize,mean.path.length=round(mean(path.length,na.rm=TRUE)),
      mean.k.edge=round(mean(K.edge,na.rm=TRUE)),
      mean.min.degree=round(mean(min_deg,na.rm=TRUE)),
      mean.max.degree=round(mean(max_deg,na.rm=TRUE)),
      mean.mean.degree=round(mean(mean_deg,na.rm=TRUE)),
      max.adjacency=max(adjacency,na.rm = TRUE)) |> 
  map(arrange,drug_pair)

# Add rates.
rates <- read.csv(file.path("data/4.PhylogeneticComparativeMethods/dataset1_ratios.csv")) 

r <- rates |>
  select(clusters,sigma.rate) |>
  distinct() |>
  arrange(sigma.rate)
clusters <- factor(r$clusters, levels = r$clusters[order(r$sigma.rate)])

df_DDI_tot <- map(DDI_dist_conn_deg_adj,inner_join,rates,by="drug_pair") |>
  map(distinct)
names(df_DDI_tot) <- net

df_DDI_tot <- bind_rows(df_DDI_tot, .id = "network") |> 
  mutate(max.adjacency=as.factor(max.adjacency)) |> 
  mutate(connection_groups=cut_interval(mean.k.edge, 6)) |> 
  mutate(connection_groups=ifelse(is.na(connection_groups),0,connection_groups)) |> 
  mutate(connection_onoff=ifelse(mean.k.edge==0,"disconnected","connected")) |> 
  mutate(int_sign_ebw = factor(int_sign_ebw,levels = c("Synergy","Additivity","Antagonism"))) |> 
  mutate(drug_cat=ifelse(drug_category=="Same",1,0)) |> 
  mutate(targ=ifelse(targeted_process=="Same",1,0)) |> 
  mutate(us=ifelse(use=="Same",1,0)) |> 
  mutate(ebw_g=ifelse(int_sign_ebw=="Antagonism",+1,
                      ifelse(int_sign_ebw=="Additivity",0,
                             ifelse(int_sign_ebw=="Synergy",-1,NA)))) |> 
  mutate(ecr_g=ifelse(int_sign_ecr=="Antagonism",+1,
                      ifelse(int_sign_ecr=="Additivity",0,
                             ifelse(int_sign_ecr=="Synergy",-1,NA)))) |> 
  mutate(seo_g=ifelse(int_sign_seo=="Antagonism",+1,
                      ifelse(int_sign_seo=="Additivity",0,
                             ifelse(int_sign_seo=="Synergy",-1,NA)))) |> 
  mutate(stm_g=ifelse(int_sign_stm=="Antagonism",+1,
                      ifelse(int_sign_stm=="Additivity",0,
                             ifelse(int_sign_stm=="Synergy",-1,NA)))) |> 
  mutate(pae_g=ifelse(int_sign_pae=="Antagonism",+1,
                      ifelse(int_sign_pae=="Additivity",0,
                             ifelse(int_sign_pae=="Synergy",-1,NA)))) |> 
  mutate(pau_g=ifelse(int_sign_pau=="Antagonism",+1,
                      ifelse(int_sign_pau=="Additivity",0,
                             ifelse(int_sign_pau=="Synergy",-1,NA)))) |> 
  rowwise() |> 
  mutate(sum_g=ebw_g+ecr_g+seo_g+stm_g+pae_g+pau_g) |> 
  distinct()


## SAVE
for (i in 1:length(DDI_dist_conn_deg_adj)){
  DDI_dist_conn_deg_adj[[i]]$network <- net[i]
  dist_conn_deg_adj[[i]]$network <- net[i]
}

df_target_tot <- map(dist_conn_deg_adj,inner_join,rates,by="drug_pair") |>
  map(distinct) |> bind_rows(.id = "network") |> 
  mutate(adjacency=as.factor(adjacency)) |> 
  mutate(connection_groups=cut_interval(K.edge, 6)) |> 
  mutate(connection_groups=ifelse(is.na(connection_groups),0,connection_groups)) |> 
  mutate(connection_onoff=ifelse(K.edge==0,"disconnected","connected"))  |> 
  mutate(int_sign_ebw = factor(int_sign_ebw,levels = c("Synergy","Additivity","Antagonism"))) |> 
  mutate(drug_cat=ifelse(drug_category=="Same",1,0)) |> 
  mutate(targ=ifelse(targeted_process=="Same",1,0)) |> 
  mutate(us=ifelse(use=="Same",1,0)) |> 
  mutate(ebw_g=ifelse(int_sign_ebw=="Antagonism",+1,
                      ifelse(int_sign_ebw=="Additivity",0,
                             ifelse(int_sign_ebw=="Synergy",-1,NA)))) |> 
  mutate(ecr_g=ifelse(int_sign_ecr=="Antagonism",+1,
                      ifelse(int_sign_ecr=="Additivity",0,
                             ifelse(int_sign_ecr=="Synergy",-1,NA)))) |> 
  mutate(seo_g=ifelse(int_sign_seo=="Antagonism",+1,
                      ifelse(int_sign_seo=="Additivity",0,
                             ifelse(int_sign_seo=="Synergy",-1,NA)))) |> 
  mutate(stm_g=ifelse(int_sign_stm=="Antagonism",+1,
                      ifelse(int_sign_stm=="Additivity",0,
                             ifelse(int_sign_stm=="Synergy",-1,NA)))) |> 
  mutate(pae_g=ifelse(int_sign_pae=="Antagonism",+1,
                      ifelse(int_sign_pae=="Additivity",0,
                             ifelse(int_sign_pae=="Synergy",-1,NA)))) |> 
  mutate(pau_g=ifelse(int_sign_pau=="Antagonism",+1,
                      ifelse(int_sign_pau=="Additivity",0,
                             ifelse(int_sign_pau=="Synergy",-1,NA)))) |> 
  rowwise() |> 
  mutate(sum_g=ebw_g+ecr_g+seo_g+stm_g+pae_g+pau_g) |>
  distinct()

df_target_tot |> colnames()
df_DDI_tot |> colnames()

write.csv(df_DDI_tot,
file.path("data/5.Targets_NetworkDistance", "df_DDI_tot_network_metrics.csv"), row.names = F)
write.csv(df_target_tot,
file.path("data/5.Targets_NetworkDistance", "df_target_tot_metrics.csv"), row.names = F)