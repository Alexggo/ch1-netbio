#!/usr/bin/Rscript

##----load---------------------------------------------------------------------
library(pacman)
p_load(igraph, tidyverse, matrixStats,future.apply,tictoc)

#Which proteins are targeted by drugs?
drug_mapping_org <- read.csv("data/3.Targets_NetworkDistance/DrugTargets3_ecoli.csv") |>
  select(Drug, KEGG_eco) |>
  distinct() |>
  arrange(Drug) |>
  filter(KEGG_eco != "")

ortho_kegg <- read.csv("data/3.Targets_NetworkDistance/ortholog_table.csv") |> 
  select(2:3,5,6,8)|> 
  mutate(KO_ID = strsplit(KO_ID, "\n")) |> 
  mutate(KEGG_eco = strsplit(eco, "\n")) |> 
  mutate(KEGG_ecr = strsplit(ecr, "\n")) |> 
  mutate(KEGG_stm = strsplit(stm, "\n")) |> 
  mutate(KEGG_pae = strsplit(pae, "\n"))  |> 
  unnest(KO_ID) |> 
  unnest(KEGG_eco) |> 
  select(KEGG_eco,KEGG_ecr,
         KEGG_stm,KEGG_pae) |> 
  mutate(KEGG_ortho=KEGG_eco)

drug_mapping_ALL <- merge(drug_mapping_org,ortho_kegg,by="KEGG_eco")

drug_mapping_eco <- drug_mapping_ALL |> 
  select(Drug,KEGG_ortho,KEGG_eco) |> 
  unnest(KEGG_eco) |> 
  mutate(KEGG=KEGG_eco) |> 
  select(Drug,KEGG,KEGG_ortho) |> 
  distinct()

drug_mapping_eco1 <- drug_mapping_eco |>
  mutate(Drug1 = Drug,
         KEGG1 = KEGG) |> 
  select(Drug1,KEGG1,KEGG_ortho) |> 
  arrange(Drug1)

drug_mapping_eco2 <- drug_mapping_eco |>
  mutate(Drug2 = Drug,
         KEGG2 = KEGG) |> 
  select(Drug2,KEGG2,KEGG_ortho)|> 
  arrange(Drug2)

drug_mapping_ecr <- drug_mapping_ALL |> 
  select(Drug,KEGG_ortho,KEGG_ecr) |> 
  unnest(KEGG_ecr)|> 
  mutate(KEGG=KEGG_ecr) |> 
  select(Drug,KEGG,KEGG_ortho) |> 
  distinct()

drug_mapping_ecr1 <- drug_mapping_ecr |>
  mutate(Drug1 = Drug,
         KEGG1 = KEGG) |> 
  select(Drug1,KEGG1,KEGG_ortho) |> 
  arrange(Drug1)

drug_mapping_ecr2 <- drug_mapping_ecr |>
  mutate(Drug2 = Drug,
         KEGG2 = KEGG) |> 
  select(Drug2,KEGG2,KEGG_ortho)|> 
  arrange(Drug2)

drug_mapping_stm <- drug_mapping_ALL |> 
  select(Drug,KEGG_ortho,KEGG_stm)|> 
  unnest(KEGG_stm)|> 
  mutate(KEGG=KEGG_stm) |> 
  select(Drug,KEGG,KEGG_ortho) |> 
  distinct()

drug_mapping_stm1 <- drug_mapping_stm |>
  mutate(Drug1 = Drug,
         KEGG1 = KEGG) |> 
  select(Drug1,KEGG1,KEGG_ortho) |> 
  arrange(Drug1)

drug_mapping_stm2 <- drug_mapping_stm |>
  mutate(Drug2 = Drug,
         KEGG2 = KEGG) |> 
  select(Drug2,KEGG2,KEGG_ortho)|> 
  arrange(Drug2)

drug_mapping_pae <- drug_mapping_ALL |>
  select(Drug,KEGG_ortho,KEGG_pae)|> 
  unnest(KEGG_pae)|> 
  mutate(KEGG=KEGG_pae) |> 
  select(Drug,KEGG,KEGG_ortho) |> 
  distinct()

drug_mapping_pae1 <- drug_mapping_pae |>
  mutate(Drug1 = Drug,
         KEGG1 = KEGG) |> 
  select(Drug1,KEGG1,KEGG_ortho) |> 
  arrange(Drug1)

drug_mapping_pae2 <- drug_mapping_pae |>
  mutate(Drug2 = Drug,
         KEGG2 = KEGG) |> 
  select(Drug2,KEGG2,KEGG_ortho)|> 
  arrange(Drug2)

d_map <- list(drug_mapping_eco,
              drug_mapping_eco,
              drug_mapping_eco,
              drug_mapping_eco,
              drug_mapping_eco,
              drug_mapping_ecr,
              drug_mapping_ecr,
              drug_mapping_pae,
              drug_mapping_pae,
              drug_mapping_stm,
              drug_mapping_stm)

all_d1_map <- list(drug_mapping_eco1,
               drug_mapping_eco1,
               drug_mapping_eco1,
               drug_mapping_eco1,
               drug_mapping_eco1,
               drug_mapping_ecr1,
               drug_mapping_ecr1,
               drug_mapping_pae1,
               drug_mapping_pae1,
               drug_mapping_stm1,
               drug_mapping_stm1)

all_d2_map <- list(drug_mapping_eco2,
                   drug_mapping_eco2,
                   drug_mapping_eco2,
                   drug_mapping_eco2,
                   drug_mapping_eco2,
                   drug_mapping_ecr2,
                   drug_mapping_ecr2,
                   drug_mapping_pae2,
                   drug_mapping_pae2,
                   drug_mapping_stm2,
                   drug_mapping_stm2)

all_targets_eco <- drug_mapping_eco$KEGG |>
  unique() |> sort()
all_targets_eco |> length()
# n!/(k!*(n-k)!)=n!/(2*n-2) possible pairwise combinations

all_targets_ecr <- drug_mapping_ecr$KEGG |>
  unique() |> sort()
all_targets_ecr |> length()

all_targets_stm <- drug_mapping_stm$KEGG |>
  unique() |> sort()
all_targets_stm |> length()

all_targets_pae <- drug_mapping_pae$KEGG |>
  unique() |> sort()
all_targets_pae |> length()

all_t <- list(all_targets_eco,
              all_targets_eco,
              all_targets_eco,
              all_targets_eco,
              all_targets_eco,
                    all_targets_ecr,
                    all_targets_ecr,
                    all_targets_pae,
                    all_targets_pae,
                    all_targets_stm,
                    all_targets_stm)

##-----------------------------------------------------------------------------
# Metabolic=1, PPI=10
#Add rates and DDI info

significant_figures <- 4

full_df_all <- read.csv(paste0("results/allddi/DDI_table_rates_allddi.csv")) |> 
  rename(clusters_all=clusters,
         sigma.rate_all=sigma.rate) 

full_df_set <- read.csv(paste0("results/sen2in1/DDI_table_rates_sen2in1.csv"))|> 
  rename(clusters_set=clusters,
         sigma.rate_set=sigma.rate) |> 
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

# Ecolinet networks
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
filepath_ecolinet <- file.path("data/3.Targets_NetworkDistance/Net_8_29_2022", filenames)
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
  lapply(select, node1, node2)

# StringDB networks.
filenamesSDB <- c("STRINGDB-E.coliK12.protein.links.detailed.v12.0.txt", "STRINGDB-E.coliK12.protein.physical.links.detailed.v12.0.txt",
                  "STRINGDB-EcoliIAI1.protein.links.detailed.v12.0.txt","STRINGDB-EcoliIAI1.protein.physical.links.detailed.v12.0.txt",
               "STRINGDB-PAO1.protein.links.detailed.v12.0.txt","STRINGDB-PAO1.protein.physical.links.detailed.v12.0.txt",
               "STRINGDB-STLT2.protein.links.detailed.v12.0.txt","STRINGDB-STLT2.protein.physical.links.detailed.v12.0.txt")

filepath_SDB <- file.path("data/3.Targets_NetworkDistance/STRINGDB", filenamesSDB)
netSDB <- c("STRINGDB-full-ebw", "STRINGDB-physical-ebw",
            "STRINGDB-full-ecr", "STRINGDB-physical-ecr",
            "STRINGDB-full-pae","STRINGDB-physical-pae",
            "STRINGDB-full-stm","STRINGDB-physical-stm")

list_netSDB <- lapply(filepath_SDB,read_table)  

# 0.9, 0.7 or 0.4 (see STRING-DB website)
list_netSDB <- list_netSDB |> 
  lapply(function(df) {
    df %>%
      filter(combined_score >= quantile(combined_score, 0.70)) %>%
      select(node1, node2)
  })

list_net <- c(list_net,list_netSDB)
net_names <- c(net,netSDB)

lapply(list_net,dim)

#Create igraphs from tables.
list_graphs <- list_net |>
  lapply(graph_from_data_frame, directed = FALSE, vertices = NULL)

list_path_conn_deg <- list()
df_join <- list()
df_join2 <- list()
for (i in seq_along(list_graphs)) {
  print(paste("current network:",i))
  # Calculate target combinations
  target_in_net <- all_t[[i]]

  T_num1 <- target_in_net |> length()
  # 27*26/2=351 possible combinations
  C_num1 <- T_num1*(T_num1-1)/2 # pairwise combinations
  C_num2 <- T_num1 # combinations with themselves
  total_c <- C_num1 + C_num2

  # next line for testing:
  #target_in_net <- target_in_net[1:2]  
  gr_vert <- t(combn(target_in_net,2))
  gr_vert_same <- cbind(target_in_net,target_in_net)
  gr_vert <- rbind(gr_vert,gr_vert_same)
  gr_vert <- gr_vert[,]
  #Expecting maximum of 351+27=378 combinations
  
  ##-----------------------------------------------------------------------------
  # Calculate path.length, k-edge connectivity,
  # and node degree for the possible targets.
  total_sample <- T_num1
  
  network <- list_graphs[[i]]
  nodes <- V(network) |> names()
  ID <- seq_along(nodes)
  # For targets
  index_targets <- cbind(nodes,ID) |> 
    as.data.frame() |> 
    mutate(ID=as.numeric(ID)) |> 
    mutate(target=ifelse(nodes%in%target_in_net,TRUE,FALSE)) |> 
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
    mutate(not_target=ifelse(nodes%in%target_in_net,FALSE,TRUE)) |> 
    arrange(desc(not_target))
  nottargets_in_net <- index_nottargets |> 
    filter(not_target==TRUE) |> 
    select(ID) |> 
    unique() |> 
    pull() |> 
    sort() |> 
    sample(size=T_num1,replace = FALSE)
  comb_nottargets <- t(combn(nottargets_in_net,2))
  comb_nottargets |> dim() #net1=45, net2=190,net3=153 target combinations
  
  comb_targets <- cbind(comb_targets,rep(1,L_targets*(L_targets-1)/2))
  comb_nottargets <- cbind(comb_nottargets,rep(0,L_targets*(L_targets-1)/2))
  comb <- rbind(comb_targets,comb_nottargets)
  
  # Calculate distance between node combinations
  print("Now doing path length")
  tic()
  allpath <- future_apply(comb,1,function(edges){
    shortest.paths(list_graphs[[i]],
                   v = edges[1], to = edges[2])
  })
  toc()

  # Calculate connectivity between node combinations
  print("Now doing connectivity")
  tic()
  kcon <- future_apply(comb,1,function(edges){
    edge_connectivity(list_graphs[[i]],
                      source = edges[1], target = edges[2])
  })
  toc()
  
  # Calculate node degree number
  unique_tar <- unique(c(comb[,1],comb[,2]))
  print("Now doing degree")
  tic()
  deg <-  igraph::degree(list_graphs[[i]],
                        v = unique_tar,mode="all")
  toc()
  
  # Calculate betweeness
  print("Now doing betweeness")
  tic()
  bet <-  igraph::betweenness(list_graphs[[i]],
                        v = unique_tar,directed=FALSE,cutoff = 15)
  toc()
  
  #Are the two nodes adjacent?
  print("Now doing adjancency")
  tic()
  adj <-  future_apply(comb,1,function(edges){
    are_adjacent(list_graphs[[i]],
                 v1 = edges[1], v2 = edges[2])
  })
  toc()
  
  # Calculate EV centrality
  print("Now doing EV centrality")
  tic()
  evcentr <- evcent(list_graphs[[i]])$vector
  toc()
  
  net_results <- data.frame(N1=comb[,1],
                            N2=comb[,2],
                            is.target=comb[,3],
                            path.length=allpath,
                            adjacency=adj,
                            K.edge=kcon)
  index_targets1 <- index_targets
  colnames(index_targets1) <- c("KEGG1","N1","is.target.1")
  index_targets2 <- index_targets
  colnames(index_targets2) <- c("KEGG2","N2","is.target.2")
  
  a <- left_join(net_results,index_targets1,by="N1")
  b <- left_join(a,index_targets2,by="N2")

  b$Between1 <- bet[b$KEGG1]
  b$Between2 <- bet[b$KEGG2]
  b$Degree1 <- deg[b$KEGG1]
  b$Degree2 <- deg[b$KEGG2]
  b$EVC1 <- evcentr[b$KEGG1]
  b$EVC2 <- evcentr[b$KEGG2]
  
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

    df_join[[i]] <- left_join(list_path_conn_deg[[i]],all_d1_map[[i]], by = "KEGG1",relationship = "many-to-many")|>
    left_join(all_d2_map[[i]], by = "KEGG2",relationship = "many-to-many") |>
    distinct() |> 
    select(-c("is.target.1","is.target.2")) |> 
    mutate(Drug1=ifelse(is.na(Drug1),"Non-targets",Drug1)) |> 
    mutate(Drug2=ifelse(is.na(Drug2),"Non-targets",Drug2))

  df_now <- df_join[[i]]
  for (ind in 1:dim(df_now)[1]){
    row_d <- df_now[ind,]
    
    current <- paste0(row_d$Drug1,"-",row_d$Drug2)
    sorted <- sort(c(row_d$Drug1,row_d$Drug2))
    sorted <- paste0(sorted[1],"-",sorted[2])
    
    if (current==sorted){
      print(paste("Row",ind,"is kept the same"))
    } else {
      A <- row_d$Drug2
      B <- row_d$Drug1
      
      df_now[ind,]$Drug1 <- A
      df_now[ind,]$Drug2 <- B
      
      A <- row_d$N2
      B <- row_d$N1
      
      df_now[ind,]$N1 <- A
      df_now[ind,]$N2 <- B
      
      A <- row_d$Degree2
      B <- row_d$Degree1
      
      df_now[ind,]$Degree1 <- A
      df_now[ind,]$Degree2 <- B
      
      A <- row_d$Between2
      B <- row_d$Between1
      
      df_now[ind,]$Between1 <- A
      df_now[ind,]$Between2 <- B
      
      A <- row_d$KEGG2
      B <- row_d$KEGG1
      
      df_now[ind,]$KEGG1 <- A
      df_now[ind,]$KEGG2 <- B
      
      A <- row_d$EVC2
      B <- row_d$EVC1
      
      df_now[ind,]$EVC1 <- A
      df_now[ind,]$EVC2 <- B
      
      print(paste("Row",ind,"is rearranged")) 
    }
  }
  df_join2[[i]] <- df_now
}

#save(df_join2,file="results/df_join2.R")
#load(file="results/df_join2.R")


df_join3 <- df_join2 |>
  lapply(mutate, DRUG_ID=paste0(Drug1,"-",Drug2)) |>
  lapply(distinct) |>
  lapply(arrange,DRUG_ID) |> 
  lapply(filter,Drug1<Drug2|Drug1=="Non-targets")

dist_conn_deg_adj <- lapply(df_join3, left_join, full_df, by = "DRUG_ID") |>
  lapply(distinct) |> 
  lapply(mutate,DRUG_ID=ifelse(DRUG_ID=="Non-targets-Non-targets","Non-targets",DRUG_ID)) |> 
  lapply(mutate,int_sign_ebw=ifelse(is.na(int_sign_ebw),"Non-targets",int_sign_ebw)) |> 
  lapply(mutate,int_sign_ecr=ifelse(is.na(int_sign_ecr),"Non-targets",int_sign_ecr)) |>
  lapply(mutate,int_sign_seo=ifelse(is.na(int_sign_seo),"Non-targets",int_sign_seo)) |>
  lapply(mutate,int_sign_stm=ifelse(is.na(int_sign_stm),"Non-targets",int_sign_stm)) |>
  lapply(mutate,int_sign_pae=ifelse(is.na(int_sign_pae),"Non-targets",int_sign_pae)) |>
  lapply(mutate,int_sign_pau=ifelse(is.na(int_sign_pau),"Non-targets",int_sign_pau)) 

names(dist_conn_deg_adj) <- net_names

df_target_tot  <- bind_rows(dist_conn_deg_adj, .id = "network") |> 
  mutate(Drug1=Drug1.x,Drug2=Drug2.x) |> 
  select(!c("Drug2.x","Drug1.x","Drug2.y","Drug1.y")) |> 
  select(network,KEGG_ID,KEGG1,KEGG2,N1,N2,is.target,
         DRUG_ID,Drug1,Drug2,
         clusters_all,sigma.rate_all,type,clusters_set,sigma.rate_set,SET,
         path.length,adjacency,K.edge,Degree1,Degree2,Between1,Between2,EVC1,EVC2,
         min_deg,max_deg,mean_deg,min_bet,max_bet,mean_bet,min_EVC,max_EVC,mean_EVC,connection_onoff,
         code_3letter1,code_3letter2,drug_cat,drug_category1,drug_category2,categorycategory,
         targ,targeted_cellular_process1,targeted_cellular_process2,processprocess,us,
         use1,use2,useuse,ebw,ecr,seo,stm,pae,pau,int_sign_ebw,int_sign_ecr,int_sign_seo,int_sign_stm,int_sign_pae,int_sign_pau,
         ebw_g,ecr_g,seo_g,stm_g,pae_g,pau_g,sum_g,
         SR_ebw_D1,SR_ecr_D1,SR_seo_D1,SR_stm_D1,SR_pae_D1,SR_pau_D1,
         SR_ebw_D2,SR_ecr_D2,SR_seo_D2,SR_stm_D2,SR_pae_D2,SR_pau_D2,
         SR_ebw,SR_ecr,SR_seo,SR_stm,SR_pae,SR_pau,
         SR_score_ebw,SR_score_ecr,SR_score_seo,SR_score_stm,SR_score_pae,SR_score_pau,SR_score_total,
         HDC_ebw_D1,HDC_ecr_D1,HDC_seo_D1,HDC_stm_D1,HDC_pae_D1,HDC_pau_D1,
         LF_ebw_D1,LF_ecr_D1,LF_seo_D1,LF_stm_D1,LF_pae_D1,LF_pau_D1,
         HDC_ebw_D2,HDC_ecr_D2,HDC_seo_D2,HDC_stm_D2,HDC_pae_D2,HDC_pau_D2,
         LF_ebw_D2,LF_ecr_D2,LF_seo_D2,LF_stm_D2,LF_pae_D2,LF_pau_D2,
         )

df_target_tot |> 
  arrange(K.edge) |> 
  arrange(DRUG_ID) |> 
  arrange(network) |>  
  write_excel_csv(paste0("results/","df_target_metrics",".csv"))

save(df_target_tot,file="results/df_target_tot.R")

# Calculate the average per drug pair.
# One target pair can have more than one drug pair.
# Similarly unique drug pair can have more than one target pair.
# Since the DDI rate data is calculated per drug combination, we can take the
# average of network metrics for each unique drug pair.

network_values_targ <- df_target_tot  |> 
  filter(DRUG_ID!="Non-targets") |> 
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
  select(network,DRUG_ID,Drug1,Drug2,
        clusters_all,sigma.rate_all,type,clusters_set,sigma.rate_set,SET,
        code_3letter1,code_3letter2,drug_cat,drug_category1,drug_category2,categorycategory,
         targ,targeted_cellular_process1,targeted_cellular_process2,processprocess,us,
         use1,use2,useuse,ebw,ecr,seo,stm,pae,pau,int_sign_ebw,int_sign_ecr,int_sign_seo,int_sign_stm,int_sign_pae,int_sign_pau,
         ebw_g,ecr_g,seo_g,stm_g,pae_g,pau_g,sum_g,
         SR_ebw_D1,SR_ecr_D1,SR_seo_D1,SR_stm_D1,SR_pae_D1,SR_pau_D1,
         SR_ebw_D2,SR_ecr_D2,SR_seo_D2,SR_stm_D2,SR_pae_D2,SR_pau_D2,
         SR_ebw,SR_ecr,SR_seo,SR_stm,SR_pae,SR_pau,
         SR_score_ebw,SR_score_ecr,SR_score_seo,SR_score_stm,SR_score_pae,SR_score_pau,SR_score_total,
         HDC_ebw_D1,HDC_ecr_D1,HDC_seo_D1,HDC_stm_D1,HDC_pae_D1,HDC_pau_D1,
         LF_ebw_D1,LF_ecr_D1,LF_seo_D1,LF_stm_D1,LF_pae_D1,LF_pau_D1,
         HDC_ebw_D2,HDC_ecr_D2,HDC_seo_D2,HDC_stm_D2,HDC_pae_D2,HDC_pau_D2,
         LF_ebw_D2,LF_ecr_D2,LF_seo_D2,LF_stm_D2,LF_pae_D2,LF_pau_D2,
         mean.path.length,mean.k.edge,mean.min.degree,mean.max.degree,
         mean.mean.degree,max.adjacency,mean.min.bet,mean.max.bet,mean.mean.bet,mean.min.EVC,
         mean.max.EVC,mean.mean.EVC,mean.sigma_all,mean.sigma_set,connection_onoff)

# df_nottarget_tot
df_nottarget_tot <- df_target_tot |> 
  filter(DRUG_ID=="Non-targets")|> 
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
  write_excel_csv(paste0("results/","df_DDI_metrics",".csv"))

save(df_DDI_tot,file="results/df_DDI_tot.R")


# plot networks
p_load(ggraph,GGally,ggnet,
       sna,network,ggplot2,patchwork,ggplotify)

adj_mat <- list_graphs |> 
  map(as_adjacency_matrix,sparse = FALSE) 

plot_list <- list()  # List to store the plots
num_targets <- list()
num_drugs <- list()
num_nodes <- list()

for (i in seq_along(list_graphs)) {
  net_graph <- network(adj_mat[[i]], directed = is_directed(list_graphs[[i]]))
  target_in_net <- all_t[[i]]
  names <- V(list_graphs[[i]]) |>names()
  print(net_names[i])
  num_targets[[i]] <- names %in% target_in_net |>
    sum()
  num_drugs[[i]] <- d_map[[i]]$Drug[target_in_net %in% names] |>
    unique() |>
    length()
  
  col_vect <- ifelse(names %in% target_in_net, "red", "lightblue")
  V(list_graphs[[i]])$color <- col_vect
  V(list_graphs[[i]])$size <- ifelse(V(list_graphs[[i]]) %in% target_in_net, 40, 5)
  V(list_graphs[[i]])$label.cex <- 0.3
  num_nodes[[i]] <- V(list_graphs[[i]]) |>length()
  mean.length <- average.path.length(list_graphs[[i]])
  
  plot_list[[i]] <- ggnet2(net_graph, label = FALSE, size = 2) +
    geom_node_point(color = 'black',
                    shape = 21, fill = col_vect, stroke = 0.5,
                    size=2)+
                    ggtitle(net_names[i])
}

pdf(paste0("results/","network_graph","_notext.pdf"),
    width=7,height=7)
plot_list
dev.off()

fig2 <- wrap_plots(plot_list[[2]],plot_list[[3]])+ 
  plot_annotation(tag_levels = 'A',tag_suffix = '.')&
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 15, hjust = 0, vjust = 2))

ggsave(
  filename = paste0("results/allddi/fig2.pdf"),
  plot = fig2,
  width = 10,
  height = 5
)

## Figures with text
plot_list <- list()  # List to store the plots
num_targets <- list()
num_drugs <- list()
num_nodes <- list()

for (i in seq_along(list_graphs)) {
  net_graph <- network(adj_mat[[i]], directed = is_directed(list_graphs[[i]]))
  target_in_net <- all_t[[i]]
  names <- V(list_graphs[[i]]) |>names()
  print(net_names[i])
  num_targets[[i]] <- names %in% target_in_net |>
    sum()
  num_drugs[[i]] <- d_map[[i]]$Drug[target_in_net %in% names] |>
    unique() |>
    length()
  
  col_vect <- ifelse(names %in% target_in_net, "red", "lightblue")
  V(list_graphs[[i]])$color <- col_vect
  V(list_graphs[[i]])$size <- ifelse(V(list_graphs[[i]]) %in% target_in_net, 40, 5)
  V(list_graphs[[i]])$label.cex <- 0.3
  num_nodes[[i]] <- V(list_graphs[[i]]) |>length()
  mean.length <- round(average.path.length(list_graphs[[i]]),3)
  
  plot_list[[i]] <- ggnet2(net_graph, label = FALSE, size = 2) +
    geom_node_point(color = 'black',
                    shape = 21, fill = col_vect, stroke = 0.5,
                    size=3)+
    annotate("text", x=0.25, y=1.2, label= net_names[[i]],size=6)+
    annotate("text", x=0.85, y=1.15, label= paste0("Average path length=",mean.length))+
    annotate("text", x=0.85, y=1.1, label= paste0("Number of nodes=",num_nodes[[i]]))+
    annotate("text", x=0.85, y=1.05, label= paste0("Number of targets=",num_targets[[i]]))+
    annotate("text", x=0.85, y=1, label= paste0("Number of drugs=",num_drugs[[i]]))+
    ggtitle(net_names[i])    
}

pdf(paste0("results/","network_graph",".pdf"),
    width=7,height=7)
plot_list
dev.off()
