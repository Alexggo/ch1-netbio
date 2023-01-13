#!/usr/bin/Rscript

##----load---------------------------------------------------------------------
library(pacman)
p_load(igraph, tidyverse, matrixStats,future.apply)


##-----------------------------------------------------------------------------
# Metabolic=1, PPI=10
#Add rates and DDI info

significant_figures <- 3

full_df_all <- read.csv(file.path("data/4.PhylogeneticComparativeMethods",
                                  "DDI_table_rates_all.csv")) |> 
  rename(clusters_all=clusters) |> 
  rename(sigma.rate_all=sigma.rate)

full_df_set <- read.csv(file.path("data/4.PhylogeneticComparativeMethods",
                                  "DDI_table_rates_set.csv"))|> 
  rename(clusters_set=clusters) |> 
  rename(sigma.rate_set=sigma.rate) |>
  select(drug_pair,clusters_set,sigma.rate_set) 

full_df <- left_join(full_df_all,full_df_set,by="drug_pair") |> 
  mutate(SET=ifelse(is.na(clusters_set),"FALSE","TRUE"))|> 
  mutate(DRUG_ID=paste0(Drug1,"-",Drug2)) |> 
  mutate(sigma.rate_all=round(sigma.rate_all,significant_figures)) |> 
  mutate(sigma.rate_set=round(sigma.rate_set,significant_figures)) |> 
  distinct() |> 
  arrange(sigma.rate_all) |> 
  mutate(clusters_all=factor(clusters_all,
                             levels = unique(full_df_all$clusters_all[order(full_df_all$sigma.rate_all)]))) |> 
  mutate(clusters_set=factor(clusters_set,
                             levels = unique(full_df_set$clusters_set[order(full_df_set$sigma.rate_set)]))) 

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

filenames <- filenames[filenames %in% 
                         c("EcoCyc.goldstandardset.csv",
                           "EcoliNet.v1.benchmark.csv",
                           "LC.INT.EcoliNet.764gene.1073link.csv")]


#Ecocyc_goldstandard:1, Small-medium_PPI:10
filepath_ecolinet <- file.path("data/5.Targets_NetworkDistance/Net_8_29_2022", filenames)
net <- c("Co-functional(EcoCyc)", "EN:allnetworks",
         "Similar genomic context", "Co-functional (GO-BP)",
         "Co-citation", "Co-expression",
         "Co-occurence of prot. domains", "Co-functional (EcoCyc/GO-BP)",
         "High-throughput PPI", "Small/medium-scale PPI",
         "Similar phylogenetic profiles")


net <- net[net %in% c("Co-functional(EcoCyc)",
                      "Co-functional (EcoCyc/GO-BP)",
                      "Small/medium-scale PPI")]

list_net <- lapply(filepath_ecolinet,read_csv) |> 
  lapply(mutate,node1=paste0("eco:", node1),
         node2=paste0("eco:", node2)) |>
  lapply(select, node1, node2)

#Create igraphs from tables.
list_graphs <- list_net |>
  lapply(graph_from_data_frame, directed = FALSE, vertices = NULL)

#Which proteins are targeted by drugs?
drug_mapping <- read.csv("data/5.Targets_NetworkDistance/DrugTargets3_ecoli.csv") |>
  select(Drug, KEGG_eco) |>
  distinct() |>
  arrange(Drug) |>
  filter(KEGG_eco != "")
all_targets <- drug_mapping$KEGG_eco |>
  unique() |> sort()
all_targets |> length()
# 27*26/2=351 possible combinations
 
#All combinations of possible target nodes
gr.vert <- t(combn(all_targets,2))
gr.vert_same <- cbind(all_targets,all_targets)
gr.vert <- rbind(gr.vert,gr.vert_same)
gr.vert <- gr.vert[,]
#Expecting maximum of 351+27=378 combinations

# plot networks
# pdf(file = file.path("networks", "network_graph.pdf"))
# 
# for (i in seq_along(list_graphs)) {
#   g1 <- list_graphs[[i]]
#   names <- V(g1) |>names()
#   print(net[i])
#   value1 <- names %in% possible_targets |>
#     sum()
#   value2 <- drug_mapping$Drug[possible_targets %in% names] |>
#     unique() |>
#     length()
#   V(g1)$color <- ifelse(names %in% possible_targets, "red", "lightblue")
#   V(g1)$size <- ifelse(V(g1) %in% possible_targets, 40, 5)
#   V(g1)$label.cex <- 0.3
#   v.number <- V(g1) |>
#     length()
#   mean.length <- average.path.length(g1)
#   plot(g1, vertex.label = "", vertex.size = 2)
#   mytitle <- paste(net[i])
#   mysubtitle1 <- paste0("Average path length = ",
#                         round(mean.length, 4), ".Number of vertexes = ", v.number)
#   mysubtitle2 <- paste0("Number of targets = ",
#                         value1, ".Number of drugs = ", value2)
#   mtext(side = 3, line = 3, at = -0.07, adj = 0, cex = 1, mytitle)
#   mtext(side = 3, line = 2, at = -0.07, adj = 0, cex = 0.7, mysubtitle1)
#   mtext(side = 3, line = 1, at = -0.07, adj = 0, cex = 0.7, mysubtitle2)
#   print(paste("Nodes", v.number, "Targets", value1, "Drugs", value2))
# }
# dev.off()
# 
# pdf(file = file.path("networks", "network_graph_notext.pdf"))
# 
# for (i in seq_along(list_graphs)) {
#   g1 <- list_graphs[[i]]
#   names <- V(g1) |>names()
#   print(net[i])
#   value1 <- names %in% possible_targets |>
#     sum()
#   value2 <- drug_mapping$Drug[possible_targets %in% names] |>
#     unique() |>
#     length()
#   V(g1)$color <- ifelse(names %in% possible_targets, "red", "lightblue")
#   V(g1)$size <- ifelse(V(g1) %in% possible_targets, 40, 5)
#   V(g1)$label.cex <- 0.3
#   v.number <- V(g1) |>length()
#   mean.length <- average.path.length(g1)
#   plot(g1, vertex.label = "", vertex.size = 2)
# }
# dev.off()

##-----------------------------------------------------------------------------
# Calculate path.length, k-edge connectivity,
# and node degree for the possible targets.
total_sample <- 4950
S_value1 <- (1+sqrt(8*total_sample+1))/2
S_value2 <- (1-sqrt(8*total_sample+1))/2
S_value <- max(S_value1,S_value2)
S_value <- round(S_value,0)
print(S_value)

list_path_conn_deg <- list()

for (i in seq_along(list_graphs)) {
  network <- list_graphs[[i]]
  nodes <- V(network) |> names()
  ID <- seq_along(nodes)
  # For targets
  index_targets <- cbind(nodes,ID) |> 
    as.data.frame() |> 
    mutate(ID=as.numeric(ID)) |> 
    mutate(target=ifelse(nodes%in%all_targets,TRUE,FALSE)) |> 
    arrange(desc(target))
  
  targets_in_net <- index_targets |> 
    filter(target==TRUE) |> 
    select(ID) |> 
    unique() |> 
    pull() |> 
    sort()
  L_targets <- length(targets_in_net)
  
  comb_targets <- t(combn(targets_in_net,2))
  comb_targets |> dim() #net1=45, net2=190,net3=153 target combinations

  # For a small sample of non-targets
  index_nottargets <- cbind(nodes,ID) |> 
    as.data.frame() |> 
    mutate(ID=as.numeric(ID)) |> 
    mutate(not_target=ifelse(nodes%in%all_targets,FALSE,TRUE)) |> 
    arrange(desc(not_target))
  nottargets_in_net <- index_nottargets |> 
    filter(not_target==TRUE) |> 
    select(ID) |> 
    unique() |> 
    pull() |> 
    sort() |> 
    sample(size=S_value,replace = FALSE)
  comb_nottargets <- t(combn(nottargets_in_net,2))
  comb_nottargets |> dim() #net1=45, net2=190,net3=153 target combinations
  
  comb_targets <- cbind(comb_targets,rep(1,L_targets*(L_targets-1)/2))
  
  
  comb_nottargets <- cbind(comb_nottargets,rep(0,S_value*(S_value-1)/2))
  comb <- rbind(comb_targets,comb_nottargets)
  

  
  # Calculate distance between node combinations
  allpath <- future_apply(comb,1,function(edges){
    shortest.paths(list_graphs[[i]],
                   v = edges[1], to = edges[2])
  })
  
  # Calculate connectivity between node combinations
  kcon <- future_apply(comb,1,function(edges){
    edge_connectivity(list_graphs[[i]],
                      source = edges[1], target = edges[2])
  })
  
  # Calculate node degree number 1
  deg1 <- future_apply(comb,1,function(edges){
    igraph::degree(list_graphs[[i]],
                   v = edges[1])
  })
  
  # Calculate node degree number 1
  deg2 <- future_apply(comb,1,function(edges){
    igraph::degree(list_graphs[[i]],
                   v = edges[2])
  })
  
  # Calculate betweeness
  bet1 <-  future_apply(comb,1,function(edges){
    igraph::betweenness(list_graphs[[i]],
                        v = edges[1],directed=FALSE,cutoff = 15)
    
  })
  
  bet2 <-  future_apply(comb,1,function(edges){
    igraph::betweenness(list_graphs[[i]],
                        v = edges[2],directed=FALSE,cutoff = 15)
    
  })
  
  #Are the two nodes adjacent?
  adj <-  future_apply(comb,1,function(edges){
    are_adjacent(list_graphs[[i]],
                 v1 = edges[1], v2 = edges[2])
  })
  
  # Calculate EV centrality
  evcentr <- evcent(list_graphs[[i]])$vector
  
  net_results <- data.frame(N1=comb[,1],
                            N2=comb[,2],
                            is.target=comb[,3],
                            path.length=allpath,
                            adjacency=adj,
                            K.edge=kcon,
                            Degree1=deg1,
                            Degree2=deg2,
                            Between1=bet1,
                            Between2=bet2)
  index_targets1 <- index_targets
  colnames(index_targets1) <- c("KEGG1","N1","is.target.1")
  index_targets2 <- index_targets
  colnames(index_targets2) <- c("KEGG2","N2","is.target.2")
  
  a <- left_join(net_results,index_targets1,by="N1")
  b <- left_join(a,index_targets2,by="N2")
  
  # Add EVC
  for (j in 1:dim(b)[1]){
    K1 <- b[j,11]
    b[j,15] <- evcentr[K1]
    
    K2 <- b[j,13]
    b[j,16] <- evcentr[K2]
  }
  colnames(b)[15] <- "EVC1"
  colnames(b)[16] <- "EVC2"
  
  list_path_conn_deg[[i]] <- b |> 
    mutate(Between1=round(Between1,2),
           Between2=round(Between2,2),
           EVC1=round(EVC1,2),
           EVC2=round(EVC2,2)) |> 
    rowwise() |> 
    mutate(min_deg=min(Degree1,Degree2,na.rm=TRUE)) |> 
    mutate(max_deg=max(Degree1,Degree2,na.rm=TRUE)) |> 
    mutate(mean_deg=(Degree1+Degree2)/2) |> 
    mutate(min_bet=min(Between1,Between2,na.rm=TRUE)) |> 
    mutate(max_bet=max(Between1,Between2,na.rm=TRUE)) |> 
    mutate(mean_bet=(Between1+Between2)/2) |>  
    mutate(min_EVC=min(EVC1,EVC2,na.rm=TRUE)) |> 
    mutate(max_EVC=max(EVC1,EVC2,na.rm=TRUE)) |> 
    mutate(mean_EVC=(EVC1+EVC2)/2) |> 
    mutate(connection_onoff =
             ifelse(K.edge == 0, "disconnected", "connected")) |> 
    ungroup() |> 
    distinct() |> 
    mutate(KEGG_ID=paste0(KEGG1,"_",KEGG2))
}

drug_mapping1 <- drug_mapping |>
  mutate(Drug1 = Drug,
         KEGG1 = KEGG_eco) |> 
  select(Drug1,KEGG1) |> 
  arrange(Drug1)

drug_mapping2 <- drug_mapping |>
  mutate(Drug2 = Drug,
         KEGG2 = KEGG_eco) |> 
  select(Drug2,KEGG2)|> 
  arrange(Drug2)

df_join <- list_path_conn_deg |>
  lapply(left_join,drug_mapping1, by = "KEGG1")

df_join1 <- df_join |>
  lapply(left_join,drug_mapping2, by = "KEGG2") |>
  lapply(distinct) |> 
  lapply(select,-c("is.target.1","is.target.2")) |> 
  lapply(mutate,Drug1=ifelse(is.na(Drug1),"Sample_nodes",Drug1)) |> 
  lapply(mutate,Drug2=ifelse(is.na(Drug2),"Sample_nodes",Drug2))

df_join2 <- list()
for (j in seq_along(df_join1)){
  
  df_now <- df_join1[[j]]
  print(paste("Doing list",j))
  
  for (i in 1:dim(df_now)[1]){
    z <- df_now[i,]
    
    current <- paste0(z$Drug1,"-",z$Drug2)
    sorted <- sort(c(z$Drug1,z$Drug2))
    sorted <- paste0(sorted[1],"-",sorted[2])
    
    if (current!=sorted){
      A <- z$Drug2
      B <- z$Drug1
      
      df_now[i,26] <- A
      df_now[i,27] <- B
      
      A <- z$N2
      B <- z$N1
      
      df_now[i,1] <- A
      df_now[i,2] <- B
      
      A <- z$Degree2
      B <- z$Degree1
      
      df_now[i,7] <- A
      df_now[i,8] <- B
      
      A <- z$Between2
      B <- z$Between1
      
      df_now[i,9] <- A
      df_now[i,10] <- B
      
      A <- z$KEGG2
      B <- z$KEGG1
      
      df_now[i,11] <- A
      df_now[i,12] <- B
      
      A <- z$EVC2
      B <- z$EVC1
      
      df_now[i,13] <- A
      df_now[i,14] <- B
      
      print(paste("Row",i,"is rearranged"))
      
    }else{
      print(paste("Row",i,"is kept the same"))
    }
  }
  df_join2[[j]] <- df_now
  print(paste("Done with list",j))
}

df_join3 <- df_join2 |>
  lapply(mutate, DRUG_ID=paste0(Drug1,"-",Drug2)) |>
  lapply(distinct) |>
  lapply(arrange,DRUG_ID) |> 
  lapply(filter,Drug1<Drug2|Drug1=="Sample_nodes")

dist_conn_deg_adj <- lapply(df_join3, left_join, full_df, by = "DRUG_ID") |>
  lapply(distinct) |> 
  lapply(mutate,DRUG_ID=ifelse(DRUG_ID=="Sample_nodes-Sample_nodes","Sample_nodes",DRUG_ID)) |> 
  lapply(mutate,int_sign_ebw=ifelse(is.na(int_sign_ebw),"Sample_nodes",int_sign_ebw)) |> 
  lapply(mutate,int_sign_ecr=ifelse(is.na(int_sign_ecr),"Sample_nodes",int_sign_ecr)) |>
  lapply(mutate,int_sign_seo=ifelse(is.na(int_sign_seo),"Sample_nodes",int_sign_seo)) |>
  lapply(mutate,int_sign_stm=ifelse(is.na(int_sign_stm),"Sample_nodes",int_sign_stm)) |>
  lapply(mutate,int_sign_pae=ifelse(is.na(int_sign_pae),"Sample_nodes",int_sign_pae)) |>
  lapply(mutate,int_sign_pau=ifelse(is.na(int_sign_pau),"Sample_nodes",int_sign_pau)) 
names(dist_conn_deg_adj) <- net

df_target_tot  <- bind_rows(dist_conn_deg_adj, .id = "network") |> 
  mutate(Drug1=Drug1.x,Drug2=Drug2.x) |> 
  select(!c("Drug2.x","Drug1.x","Drug2.y","Drug1.y")) |> 
  select(network,KEGG_ID,KEGG1,KEGG2,N1,N2,is.target,
         DRUG_ID,Drug1,Drug2,112:117,
         5:11,14:25,30:111)

df_target_tot |> 
  arrange(K.edge) |> 
  arrange(DRUG_ID) |> 
  arrange(network) |>  
  write_excel_csv2(
    file.path("data/5.Targets_NetworkDistance",
              paste0("df_target_metrics.csv")))

# Calculate the average per drug pair.
# One target pair can have more than one drug pair.
# Similarly unique drug pair can have more than one target pair.
# Since the DDI rate data is calculated per drug combination, we can take the
# average of network metrics for each unique drug pair.

network_values_targ <- df_target_tot  |> 
  filter(DRUG_ID!="Sample_nodes") |> 
  group_by(network,DRUG_ID) |> 
  summarise(mean.path.length = ifelse(mean(path.length,na.rm=TRUE)==Inf,Inf,round(mean(path.length,na.rm=TRUE),significant_figures)),
            mean.k.edge = round(mean(K.edge, na.rm = TRUE),significant_figures),
            mean.min.degree = round(mean(min_deg, na.rm = TRUE),significant_figures),
            mean.max.degree = round(mean(max_deg, na.rm = TRUE),significant_figures),
            mean.mean.degree = round(mean(mean_deg, na.rm = TRUE),significant_figures),
            max.adjacency = max(adjacency, na.rm = TRUE),
            mean.min.bet = round(mean(min_bet, na.rm = TRUE),significant_figures),
            mean.max.bet = round(mean(max_bet, na.rm = TRUE),significant_figures),
            mean.mean.bet = round(mean(mean_bet, na.rm = TRUE),significant_figures),
            mean.min.EVC = round(mean(min_EVC, na.rm = TRUE),significant_figures),
            mean.max.EVC = round(mean(max_EVC, na.rm = TRUE),significant_figures),
            mean.mean.EVC = round(mean(mean_EVC, na.rm = TRUE),significant_figures),
            mean.sigma_all = round(mean(sigma.rate_all, na.rm = TRUE),significant_figures),
            mean.sigma_set = round(mean(sigma.rate_set, na.rm = TRUE),significant_figures)) |> 
  ungroup() |> 
  distinct()

df_DDI_targ  <- right_join(full_df,network_values_targ, by = "DRUG_ID") |>
  distinct() |>
  mutate(max.adjacency = as.factor(max.adjacency)) |>
  mutate(connection_onoff =
           ifelse(mean.k.edge == 0, "disconnected", "connected")) |> 
  ungroup() |> 
  distinct() |> 
  filter(!is.na(Drug1))|> 
  arrange(mean.k.edge) |> 
  arrange(DRUG_ID) |> 
  select(network,DRUG_ID,Drug1,Drug2,SET,
         87:91,5:86,95:109)

# df_nottarget_tot
df_nottarget_tot <- df_target_tot |> 
  filter(DRUG_ID=="Sample_nodes")|> 
  mutate(DRUG_ID=KEGG_ID)

network_values_nottarg <- df_nottarget_tot |> 
  group_by(DRUG_ID) |> 
  summarise(mean.path.length = ifelse(mean(path.length,na.rm=TRUE)==Inf,Inf,round(mean(path.length,na.rm=TRUE),significant_figures)),
            mean.k.edge = round(mean(K.edge, na.rm = TRUE),significant_figures),
            mean.min.degree = round(mean(min_deg, na.rm = TRUE),significant_figures),
            mean.max.degree = round(mean(max_deg, na.rm = TRUE),significant_figures),
            mean.mean.degree = round(mean(mean_deg, na.rm = TRUE),significant_figures),
            max.adjacency = max(adjacency, na.rm = TRUE),
            mean.min.bet = round(mean(min_bet, na.rm = TRUE),significant_figures),
            mean.max.bet = round(mean(max_bet, na.rm = TRUE),significant_figures),
            mean.mean.bet = round(mean(mean_bet, na.rm = TRUE),significant_figures),
            mean.min.EVC = round(mean(min_EVC, na.rm = TRUE),significant_figures),
            mean.max.EVC = round(mean(max_EVC, na.rm = TRUE),significant_figures),
            mean.mean.EVC = round(mean(mean_EVC, na.rm = TRUE),significant_figures),
            mean.sigma_all = round(mean(sigma.rate_all, na.rm = TRUE),significant_figures),
            mean.sigma_set = round(mean(sigma.rate_set, na.rm = TRUE),significant_figures)) |> 
  ungroup() |> 
  distinct() 


df_DDI_nottarg  <- right_join(df_nottarget_tot,network_values_nottarg, by = "DRUG_ID") |>
  distinct() |>
  mutate(max.adjacency = as.factor(max.adjacency)) |>
  mutate(connection_onoff =
           ifelse(mean.k.edge == 0, "disconnected", "connected")) |> 
  ungroup() |> 
  distinct() |> 
  arrange(mean.k.edge) |> 
  arrange(DRUG_ID) |> 
  select(colnames(df_DDI_targ))

#rbind df_DDI_targ and df_DDI_nottarg
df_DDI_tot <- rbind(df_DDI_targ,df_DDI_nottarg)


df_DDI_tot    |> 
  write_excel_csv2(
    file.path("data/5.Targets_NetworkDistance",
              paste0("df_DDI_metrics.csv")))


