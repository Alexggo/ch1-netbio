#' ---
#' title: "5.0.Standardize networks"
#' author: "Alex Gil"
#' 
library(pacman)
p_load(igraph,tidyverse,BioNet,ape,geomorph)


# Convert enzymes to kegg ids for metabolic network.
edges <- read_delim("data/5.Targets_NetworkDistance/MODEL1108160000_BASE.txt",col_names = TRUE)
node_info <- read_csv("data/5.Targets_NetworkDistance/MODEL1108160000_BASE_node.csv") %>% 
  select(name,GENE_ASSOCIATION)

edges %>% head()
node_info %>% head()

node_info <- node_info %>%  mutate(GENE=gsub("\\(","",as.character(GENE_ASSOCIATION))) %>% 
  mutate(GENE=gsub("\\)","",as.character(GENE))) %>% 
  mutate(GENE=gsub("and","or",as.character(GENE))) %>% 
  separate(GENE,sep="or",into=c("c1","c2","c3","c4",
                                "c5","c6","c7","c8",
                                "c9","c10","c11","c12",
                                "c13","c14","c15","c16",
                                "c17","c18")) %>% 
  select(-GENE_ASSOCIATION) %>% 
  pivot_longer(names_to = "col",values_to = "gene",2:19) %>% 
  drop_na() %>% 
  select(name,gene) %>% 
  mutate(gene=gsub(" ","",as.character(gene))) %>% 
  distinct() %>% 
  rename(node1=name,gene1=gene) %>% 
  mutate(gene1=paste0("eco:",gene1))

d1 <- left_join(edges,node_info,by="node1")
node_info <- node_info %>% rename(node2=node1,
                                  gene2=gene1)
d2 <- left_join(d1,node_info,by="node2")
d3 <- d2 %>% 
  mutate(gene1=ifelse(is.na(gene1),node1,gene1)) %>% 
  mutate(gene2=ifelse(is.na(gene2),node2,gene2)) %>% 
  select(gene1,gene2) %>% 
  rename(node1=gene1,node2=gene2)
write_csv(d3,"data/5.Targets_NetworkDistance/MODEL1108160000_edgelist.csv")

# Fix metabolic network 2.
edges <- read_delim("data/5.Targets_NetworkDistance/MODEL1302140002_BASE.txt")
node_info <- read_csv("data/5.Targets_NetworkDistance/MODEL1302140002_BASE_node.csv") %>% 
  select(name,kegg.genes,kegg.compound) %>% 
  rowwise() %>% 
  mutate(kegg=ifelse(!is.na(kegg.genes),kegg.genes,name)) %>% 
  select(name,kegg) %>% 
  distinct() %>% 
  rename(node1=name,gene1=kegg)

edges %>% head()
node_info %>% head()

d1 <- left_join(edges,node_info,by="node1")
node_info <- node_info %>% rename(node2=node1,
                                  gene2=gene1)
d2 <- left_join(d1,node_info,by="node2")
d3 <- d2 %>% 
  mutate(gene1=ifelse(is.na(gene1),node1,gene1)) %>% 
  mutate(gene2=ifelse(is.na(gene2),node2,gene2)) %>% 
  select(gene1,gene2) %>% 
  rename(node1=gene1,node2=gene2)
write_csv(d3,"data/5.Targets_NetworkDistance/MODEL130214000_edgelist.csv")


#Standardize networks to the same format.

filename_ecolinet <- c("EcoCyc.goldstandardset.txt","EcoliNet.v1.txt","GN.INT.EcoliNet.3568gene.23439link.txt",
              "GO-BP.goldstandardset.txt","CC.EcoliNet.v1.2296gene.50528link.txt","CX.INT.EcoliNet.v1.4039gene.67494link.txt",
              "DC.EcoliNet.2283gene.9643link.txt","EcoliNet.v1.benchmark.txt","HT.INT.EcoliNet.3209gene.15543link.txt",
              "LC.INT.EcoliNet.764gene.1073link.txt","PG.INT.EcoliNet.v1.1817gene.17504link.txt")

filename_met <- c("MODEL1108160000_edgelist.csv","MODEL130214000_edgelist.csv")

filenames <- c(filename_ecolinet,filename_met)

list_ecolinet <- map(file.path("data/5.Targets_NetworkDistance",filename_ecolinet),read_table) 
list_met <- map(file.path("data/5.Targets_NetworkDistance",filename_met),read_csv)
# Add eco: to ecolinet networks.
for (i in 1:length(list_ecolinet)){
  list_ecolinet[[i]]$node1 <- paste0("eco:",list_ecolinet[[i]]$node1)
  list_ecolinet[[i]]$node2 <- paste0("eco:",list_ecolinet[[i]]$node2)
}
list_net <- c(list_ecolinet,list_met)

# Which proteins are targeted by drugs?
nodes <- read.csv("data/5.Targets_NetworkDistance/DrugTargets3_ecoli.csv") %>% 
  select(Drug,KEGG_eco) %>% distinct() %>% arrange(Drug) %>% filter(KEGG_eco!="")
possible_targets <- nodes$KEGG_eco %>% unique()
# Calculate path.length, k-edge and node degree for these targets.
list_net <- list_net %>% map(select,node1,node2) %>% 
  map(rename,Uni_N1=node1) %>% 
  map(rename,Uni_N2=node2)
list_graphs <- list_net %>% 
  map(graph_from_data_frame,directed=FALSE,vertices=NULL)
l_deg <- list_graphs %>% 
  map(igraph::degree)
l_deg <- list_graphs %>% 
  map(igraph::degree)
list_degdist <- list_graphs %>% 
  map(degree_distribution) %>% 
  map(plot)
mean_pathlength <- list_graphs %>% 
  map_dbl(average.path.length)
shortest_paths <-  list_graphs %>% 
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
  drugs <- nodes %>% filter(KEGG_eco%in%targets) %>% select(Drug) %>% unique() %>% pull() %>% sort()
  drug_list[[i]] <- drugs
  print(paste("Network:",filenames[[i]]))
  Tar <-   rownames(dist_mat) %in% possible_targets %>% sum() 
  print(paste("Number of targets:",Tar))
  print(targets)
  D <-   drugs %>% length()
  print(paste("Number of drugs:",D))
  drugs %>% print()
  # Small matrix of targets from Brochado
  distnew2 <- dist_mat[r,c] %>% 
    as.data.frame()
  # Convert matrix to table
  distnew2 <- distnew2[ , order(names(distnew2))]
  distnew2$KEGG1 <- rownames(distnew2)
  distnew2 <- distnew2 %>%
    arrange(KEGG1) %>% 
    select(-KEGG1)
  distnew2[lower.tri(distnew2)] <- NA
  distnew2$KEGG1 <- colnames(distnew2)
  distnew2 <- distnew2 %>%
    select(dim(distnew2)[2],1:(dim(distnew2)[2]-1))
  list_path[[i]] <- distnew2  %>% 
    pivot_longer(names_to="KEGG2",values_to="path.length",2:(dim(distnew2)[2]))%>%
    filter(!is.na(path.length))
}


#k.edge.connectivity between known nodes that are in g1
list_Kedge <- list()
for (i in 1:length(list_graphs)){
  selection <- possible_targets[possible_targets %in% rownames(shortest_paths[[i]])]
  t <- matrix(nrow = length(selection),ncol=length(selection))
  for (j in 1:length(selection)){
    for (k in 1:length(selection)){
      if(j==k){t[j,k] <- NA}else{
        t[j,k] <- edge_connectivity(list_graphs[[i]],source = selection[j],target = selection[k])
        colnames(t) <- selection
        rownames(t) <- selection
        
        connnew <- t %>% as.data.frame()
        connnew <- connnew[ , order(names(connnew))]
        connnew$KEGG1 <- rownames(connnew)
        connnew <- connnew %>%
          arrange(KEGG1) %>% 
          select(-KEGG1)
        connnew[lower.tri(connnew)] <- NA
        connnew$KEGG1 <- colnames(connnew)
        connnew1 <- connnew %>%
          select(dim(connnew)[2],1:(dim(connnew)[2]-1))
        connnew2 <- connnew1  %>% 
          pivot_longer(names_to="KEGG2",values_to="K.edge",2:(dim(connnew)[2]))%>%
          filter(!is.na(K.edge)) %>% 
          filter(K.edge!=0) 
        list_Kedge[[i]] <- connnew2
      }
    }
  }
}

# Merge both path length and K edge connectivity.
list_path <- list_path %>% 
  map(mutate,KEGG1_KEGG2=paste0(KEGG1,"_",KEGG2))
list_Kedge <- list_Kedge %>% 
  map(mutate,KEGG1_KEGG2=paste0(KEGG1,"_",KEGG2))

dist_conn <- map2(list_path, list_Kedge, full_join, by = "KEGG1_KEGG2") %>% 
  map(select,KEGG1.x,KEGG2.x,path.length,K.edge) %>% 
  map(rename,KEGG1=KEGG1.x,
      KEGG2=KEGG2.x,
      path.length=path.length) %>% 
  map(arrange,K.edge) %>% 
  map(distinct) %>% 
  map(filter,KEGG1!=KEGG2)
 

l_deg_df <- l_deg %>% map(as.data.frame)

for (i in 1:length(l_deg_df)){
  KEGG <- l_deg_df[[i]] %>% rownames()
  degree <- l_deg_df[[i]][,1]
  l_deg_df[[i]] <- data.frame("KEGG1"=KEGG,
                              "Degree1"=degree) %>% 
    distinct()
}

tabnew2 <- map2(dist_conn, l_deg_df, left_join, by = "KEGG1") 
l_deg_df2 <- l_deg_df
l_deg_df2 <- l_deg_df2 %>% map(rename,"KEGG2"=KEGG1,"Degree2"=Degree1)
tabnew2 <- map2(tabnew2, l_deg_df2, left_join, by = "KEGG2") %>% 
  map(mutate,KEGG1_KEGG2=paste0(KEGG1,"_",KEGG2))

sum_tab2 <- tabnew2 %>% 
  map(pivot_longer,names_to="category",values_to="degree",5:6) %>% 
  map(group_by,KEGG1,KEGG2) %>% 
  map(summarise,mean_deg=mean(degree),#maybe round means.
      min_deg=min(degree),
      max_deg=max(degree))  %>% 
  map(mutate,KEGG1_KEGG2=paste0(KEGG1,"_",KEGG2)) %>% 
  map(ungroup) %>% 
  map(distinct)

tabnew2 <- map2(tabnew2, sum_tab2, left_join, by = "KEGG1_KEGG2") %>% 
  map(select,-c(KEGG1.y,KEGG2.y)) %>% 
  map(rename,KEGG1=KEGG1.x,KEGG2=KEGG2.x)

nodes1 <- nodes %>%
  rename(Drug1=Drug,
         KEGG1=KEGG_eco)

df_join <- tabnew2 %>% 
  map(full_join,nodes1,by="KEGG1")

nodes2 <- nodes %>%
  rename(Drug2=Drug,
         KEGG2=KEGG_eco)
df_join1 <- df_join %>% 
  map(full_join,nodes2,by="KEGG2") %>% 
  map(mutate,KEGG1_KEGG2=paste0(KEGG1,"_",KEGG2)) %>% 
  map(distinct) %>% 
  map(mutate,Drug1_Drug2=paste0(Drug1,"_",Drug2)) %>% 
  map(select,KEGG1_KEGG2,KEGG1,KEGG2,Drug1_Drug2,Drug1,Drug2,path.length,K.edge,Degree1,Degree2,
      mean_deg.x,min_deg.x,max_deg.x) %>% 
  map(rename,mean_deg=mean_deg.x,min_deg=min_deg.x,max_deg=max_deg.x) %>% 
  map(drop_na) %>% 
  map(distinct) %>% 
  map(filter,Drug1!=Drug2) %>% 
  map(filter,Drug1<Drug2)

# A DDI can have several target_target combinations. We take the average of all metrics per DDI.
df_join_sum <- df_join1 %>% 
  map(group_by,Drug1_Drug2) %>% 
  map(summarise,mean.path.length=round(mean(path.length,na.rm=TRUE)),
                mean.k.edge=round(mean(K.edge,na.rm=TRUE)),
                mean.min.degree=round(mean(min_deg,na.rm=TRUE)),
                mean.max.degree=round(mean(max_deg,na.rm=TRUE)),
                mean.mean.degree=round(mean(mean_deg,na.rm=TRUE))) %>% 
  map(distinct) %>% 
  map(arrange,Drug1_Drug2)

for (i in 1:length(df_join_sum)){
  df_join1[[i]]$network <- filenames[[i]]
  df_join_sum[[i]]$network <- filenames[[i]]
}






rates <- read.csv(file.path("data/4.PhylogeneticComparativeMethods/dataset1_ratios.csv")) %>% 
  rename(Drug1_Drug2=drug_pair) %>% 
  filter(drug1>drug2)



df_join1 <- do.call("rbind", df_join1)
df_join_sum <- do.call("rbind", df_join_sum)

write.csv(df_join1,file.path("data/5.Targets_NetworkDistance","netvalues.csv"),row.names = F)
write.csv(df_join_sum,file.path("data/5.Targets_NetworkDistance","netvalues_target.csv"),row.names = F)


pdf(file = file.path("data/5.Targets_NetworkDistance","network_graphs.pdf"))

for (i in 1:length(list_graphs)){
  g1 <- list_graphs[[i]]
  names <- V(g1) %>% names()
  print(filenames[i])
  value1 <- names %in% nodes$KEGG_eco %>% sum()
  value2 <- nodes$Drug[nodes$KEGG_eco%in%names] %>% unique() %>% length()
  V(g1)$color <- ifelse(names %in% nodes$KEGG_eco, "red", "lightblue")
  V(g1)$size <- ifelse(V(g1) %in% nodes$KEGG_eco, 40, 5)
  V(g1)$label.cex <- 0.3
  v.number <- V(g1) %>% length()
  mean.length <- average.path.length(g1)
  plot(g1,vertex.label="",vertex.size=2)
  mytitle = paste(filenames[i])
  mysubtitle1 = paste0("Average path length =",round(mean.length,4),". Number of vertexes =",v.number)
  mysubtitle2 = paste0("Number of targets =",value1,". Number of drugs =",value2)
  mtext(side=3, line=3, at=-0.07, adj=0, cex=1, mytitle)
  mtext(side=3, line=2, at=-0.07, adj=0, cex=0.7, mysubtitle1)
  mtext(side=3, line=1, at=-0.07, adj=0, cex=0.7, mysubtitle2)
}
dev.off()

pdf(file = file.path("data/5.Targets_NetworkDistance","network_graphs_notext.pdf"))

for (i in 1:length(list_graphs)){
  g1 <- list_graphs[[i]]
  names <- V(g1) %>% names()
  print(filenames[i])
  value1 <- names %in% nodes$KEGG_eco %>% sum()
  value2 <- nodes$Drug[nodes$KEGG_eco%in%names] %>% unique() %>% length()
  V(g1)$color <- ifelse(names %in% nodes$KEGG_eco, "red", "lightblue")
  V(g1)$size <- ifelse(V(g1) %in% nodes$KEGG_eco, 40, 5)
  V(g1)$label.cex <- 0.3
  v.number <- V(g1) %>% length()
  mean.length <- average.path.length(g1)
  plot(g1,vertex.label="",vertex.size=2)
}
dev.off()
