#' ---
#' title: "5.1.Calculate network parameters"
#' author: "Alex Gil"
#' 
library(igraph)
library(tidyverse)
library(BioNet)
library(ape)
library(geomorph)

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
#- MODEL1108160000_edgelist.txtnet.txt

filename <- c("EcoCyc.goldstandardset.txt","EcoliNet.v1.txt","GN.INT.EcoliNet.3568gene.23439link.txt",
              "GO-BP.goldstandardset.txt","CC.EcoliNet.v1.2296gene.50528link.txt","CX.INT.EcoliNet.v1.4039gene.67494link.txt",
              "DC.EcoliNet.2283gene.9643link.txt","EcoliNet.v1.benchmark.txt","HT.INT.EcoliNet.3209gene.15543link.txt",
              "LC.INT.EcoliNet.764gene.1073link.txt","PG.INT.EcoliNet.v1.1817gene.17504link.txt"
              ,"MODEL1108160000_edgelist.txt"
)


# Map drugs to targets (for those that are known)
nodes <- read.csv("data/5.Targets_NetworkDistance/DrugTargets3_ecoli.csv")
nodes1 <- nodes %>%
  select(KEGG_eco)%>%
  unique() %>%
  pull()

list.graph <- list()
list.results <- list()


for (i in 1:length(filename)){
  t <- file.path("data/5.Targets_NetworkDistance",paste0(filename[i],"net.txt"))
  x1 <- read_csv(t)
  x1 <- x1 %>% select(Uni_N1,Uni_N2)
  g1 <- graph_from_data_frame(x1, directed = FALSE, vertices = NULL)
  list.graph[[i]] <- g1
  deg <-  igraph::degree(g1)
  degree_distribution(g1) %>% plot()
  mean.length <- average.path.length(g1)
  #Calculate shortest path
  dist_mat <- shortest.paths(g1)
  print(filename[i])
  # Filter columns and rows for known targets
  r <- rownames(dist_mat) %in% nodes1
  c <- colnames(dist_mat) %in% nodes1
  # How many targets have drugs in Brochado?
  print("Number of targets in Brochado")
  rownames(dist_mat) %in% nodes1 %>% sum()
  # Small matrix of targets from Brochado
  distnew2 <- dist_mat[r,c] %>% 
    as.data.frame()
  # Convert matrix to table
  distnew2 <- distnew2[ , order(names(distnew2))]
  distnew2$Uniprot1 <- rownames(distnew2)
  distnew2 <- distnew2 %>%
    arrange(Uniprot1) %>% 
    select(-Uniprot1)
  distnew2[lower.tri(distnew2)] <- NA
  distnew2$Uniprot1 <- colnames(distnew2)
  distnew2 <- distnew2 %>%
    select(dim(distnew2)[2],1:(dim(distnew2)[2]-1))
  distnew2 <- distnew2  %>% 
    pivot_longer(names_to="Uniprot2",values_to="path.length",2:(dim(distnew2)[2]))%>%
    filter(!is.na(path.length))
  
  #k.edge.connectivity for known nodes that are in g1
  selection <- nodes1[nodes1 %in% rownames(dist_mat)]
  t <- matrix(nrow = length(selection),ncol=length(selection))
  for (j in 1:length(selection)){
    for (k in 1:length(selection)){
      if(j==k){t[j,k] <- NA}else{
        t[j,k] <- edge_connectivity(g1,source = selection[j],target = selection[k])
      }
    }
  }
  
  colnames(t) <- selection
  rownames(t) <- selection
  
  # Convert K matrix to table
  connnew <- t %>% as.data.frame()
  connnew <- connnew[ , order(names(connnew))]
  connnew$Uniprot1 <- rownames(connnew)
  connnew <- connnew %>%
    arrange(Uniprot1) %>% 
    select(-Uniprot1)
  connnew[lower.tri(connnew)] <- NA
  connnew$Uniprot1 <- colnames(connnew)
  connnew1 <- connnew %>%
    select(dim(connnew)[2],1:(dim(connnew)[2]-1))
  connnew2 <- connnew1  %>% 
    pivot_longer(names_to="Uniprot2",values_to="K.edge",2:(dim(connnew)[2]))%>%
    filter(!is.na(K.edge)) %>% 
    filter(K.edge!=0) 
  # Merge distance and connectivity results
  distnew2 <- distnew2 %>% 
    mutate(U_ID=paste0(Uniprot1,"-",Uniprot2))
  connnew2 <- connnew2 %>% 
    mutate(U_ID=paste0(Uniprot1,"-",Uniprot2))
  
  dist_conn <- full_join(distnew2,connnew2,by="U_ID")
  dist_conn <- dist_conn %>% 
    select(Uniprot1.x,Uniprot2.x,path.length,K.edge) %>% 
    rename(Uniprot1=Uniprot1.x,
           Uniprot2=Uniprot2.x,
           path.length=path.length) %>% 
    arrange(K.edge)
  
  # Add the degree for each node, and calculate the min, mean and max for each pair.
  tab_deg1 <- data.frame("Uniprot1"=names(deg),
                         "Degree1"=deg)
  rownames(tab_deg1) <- NULL
  
  tabnew2 <- inner_join(dist_conn,tab_deg1,by="Uniprot1")
  tab_deg2 <- tab_deg1
  colnames(tab_deg2) <- c("Uniprot2","Degree2")
  tabnew2<- inner_join(tabnew2,tab_deg2,by="Uniprot2")
  
  tabnew2$mean_deg <- apply(tabnew2[,5:6],1,mean)
  tabnew2$min_deg <- apply(tabnew2[,5:6],1,min)
  tabnew2$max_deg <- apply(tabnew2[,5:6],1,max)
  tabnew2 <- tabnew2 %>% mutate(mean_deg=round(mean_deg)) %>% 
    distinct()
  
  # Map Targets to Drugs. 
  # One target can have multiple drugs.
  unique_targets <- c(tabnew2$Uniprot1,tabnew2$Uniprot2) %>%
    unique()
  drug_target <-nodes %>%
    filter(KEGG_eco %in% unique_targets)
  
  drug_target1 <- drug_target %>% 
    select(Drug,KEGG_eco)
  colnames(drug_target1) <- c("drug1","Uniprot1")
  print("How many drugs?")
  drug_target1$drug1 %>% unique() %>% sort()
  
  tab1 <- inner_join(tabnew2,drug_target1,by="Uniprot1") %>% 
    distinct()
  colnames(drug_target1) <- c("drug2","Uniprot2")
  tab1 <- inner_join(tab1,drug_target1,by="Uniprot2")
  
  tab1 <- tab1 %>% filter(tab1$drug1!=tab1$drug2)%>% distinct()
  
  for (l in 1:dim(tab1)[1]){
    nam <- sort(c(tab1$drug1[l],tab1$drug2[l]),decreasing = F)
    
    tab1$drug_pair[l] <- paste0(nam, collapse="_")
  }
  tab1 <- tab1 %>% select(-c(drug1,drug2)) %>% distinct()
  y <- tab1 %>% group_by(drug_pair) %>% 
    summarise(mean.path.length=round(mean(path.length,na.rm=TRUE)),
              mean.k.edge=round(mean(K.edge,na.rm=TRUE)),
              mean.min.degree=round(mean(min_deg,na.rm=TRUE)),
              mean.max.degree=round(mean(max_deg,na.rm=TRUE)),
              mean.mean.degree=round(mean(mean_deg,na.rm=TRUE))) %>% 
    distinct() %>% arrange(mean.path.length) 
  y$network <- filename[i]
  tab1$network <- filename[i]
  
  list.results[[i]] <- y %>% as.data.frame()
  write.csv(y,file.path("data/5.Targets_NetworkDistance",paste0(filename[i],"netvalues.csv")),row.names = F)
  write.csv(tab1,file.path("data/5.Targets_NetworkDistance",paste0(filename[i],"netvalues_target.csv")),row.names = F)
}

pdf(file = file.path("data/5.Targets_NetworkDistance","network_graphs.pdf"))

for (i in 1:12){
  g1 <- list.graph[[i]]
  names <- V(g1) %>% names()
  print(filename[i])
  value1 <- names %in% nodes1 %>% sum()
  value2 <- nodes$Drug[nodes$KEGG_eco%in%names] %>% unique() %>% length()
  V(g1)$color <- ifelse(names %in% nodes1, "red", "lightblue")
  V(g1)$size <- ifelse(V(g1) %in% nodes1, 40, 5)
  V(g1)$label.cex <- 0.3
  v.number <- V(g1) %>% length()
  mean.length <- average.path.length(g1)
  plot(g1,vertex.label="",vertex.size=2)
  mytitle = paste(filename[i])
  mysubtitle1 = paste0("Average path length =",round(mean.length,4),". Number of vertexes =",v.number)
  mysubtitle2 = paste0("Number of targets =",value1,". Number of drugs =",value2)
  mtext(side=3, line=3, at=-0.07, adj=0, cex=1, mytitle)
  mtext(side=3, line=2, at=-0.07, adj=0, cex=0.7, mysubtitle1)
  mtext(side=3, line=1, at=-0.07, adj=0, cex=0.7, mysubtitle2)
}
dev.off()

pdf(file = file.path("data/5.Targets_NetworkDistance","network_graphs_notext.pdf"))

for (i in 1:12){
  g1 <- list.graph[[i]]
  names <- V(g1) %>% names()
  print(filename[i])
  value1 <- names %in% nodes1 %>% sum()
  value2 <- nodes$Drug[nodes$KEGG_eco%in%names] %>% unique() %>% length()
  V(g1)$color <- ifelse(names %in% nodes1, "red", "lightblue")
  V(g1)$size <- ifelse(V(g1) %in% nodes1, 40, 5)
  V(g1)$label.cex <- 0.3
  v.number <- V(g1) %>% length()
  mean.length <- average.path.length(g1)
  plot(g1,vertex.label="",vertex.size=2)
}
dev.off()
