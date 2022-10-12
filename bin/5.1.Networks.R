##----load---------------------------------------------------------------------
library(pacman)
p_load(igraph, tidyverse, ape, matrixStats,
ggpubr, plotrix, tidymodels, here,future.apply)


##-----------------------------------------------------------------------------
# Metabolic=1, PPI=10

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

#Ecocyc_goldstandard:1, Small-medium_PPI:10
filepath_ecolinet <- file.path("data/5.Targets_NetworkDistance/EcoliNet_8_29_2022", filenames)
net <- c("Co-functional(EcoCyc)", "EN:allnetworks",
"Similar genomic context", "Co-functional (GO-BP)",
"Co-citation", "Co-expression",
"Co-occurence of prot. domains", "Co-functional (EcoCyc/GO-BP)",
"High-throughput PPI", "Small/medium-scale PPI",
"Similar phylogenetic profiles")

print(paste("Filename:",filenames))
print(paste("Filepath:",filepath_ecolinet))


list_net <- future_lapply(filepath_ecolinet,read_csv) |> 
  lapply(mutate,node1=paste0("eco:", node1),
      node2=paste0("eco:", node2)) |>
  lapply(select, node1, node2)

names(list_net) <- net


#Create igraphs from tables.
list_graphs <- list_net |>
  lapply(graph_from_data_frame, directed = FALSE, vertices = NULL)

#Which proteins are targeted by drugs?
drug_mapping <- read.csv("data/5.Targets_NetworkDistance/DrugTargets3_ecoli.csv") |>
  select(Drug, KEGG_eco) |>
  distinct() |>
  arrange(Drug) |>
  filter(KEGG_eco != "")
possible_targets <- drug_mapping$KEGG_eco |>
  unique()

##-----------------------------------------------------------------------------
pdf(file = file.path("networks", "network_graph.pdf"))

for (i in seq_along(list_graphs)) {
g1 <- list_graphs[[i]]
names <- V(g1) |>names()
print(net[i])
value1 <- names %in% possible_targets |>
  sum()
value2 <- drug_mapping$Drug[possible_targets %in% names] |>
  unique() |>
  length()
V(g1)$color <- ifelse(names %in% possible_targets, "red", "lightblue")
V(g1)$size <- ifelse(V(g1) %in% possible_targets, 40, 5)
V(g1)$label.cex <- 0.3
v.number <- V(g1) |>
  length()
mean.length <- average.path.length(g1)
plot(g1, vertex.label = "", vertex.size = 2)
mytitle <- paste(net[i])
mysubtitle1 <- paste0("Average path length = ",
round(mean.length, 4), ".Number of vertexes = ", v.number)
mysubtitle2 <- paste0("Number of targets = ",
value1, ".Number of drugs = ", value2)
mtext(side = 3, line = 3, at = -0.07, adj = 0, cex = 1, mytitle)
mtext(side = 3, line = 2, at = -0.07, adj = 0, cex = 0.7, mysubtitle1)
mtext(side = 3, line = 1, at = -0.07, adj = 0, cex = 0.7, mysubtitle2)
print(paste("Nodes", v.number, "Targets", value1, "Drugs", value2))
}
dev.off()

pdf(file = file.path("networks", "network_graph_notext.pdf"))

for (i in seq_along(list_graphs)) {
g1 <- list_graphs[[i]]
names <- V(g1) |>names()
print(net[i])
value1 <- names %in% possible_targets |>
  sum()
value2 <- drug_mapping$Drug[possible_targets %in% names] |>
  unique() |>
  length()
V(g1)$color <- ifelse(names %in% possible_targets, "red", "lightblue")
V(g1)$size <- ifelse(V(g1) %in% possible_targets, 40, 5)
V(g1)$label.cex <- 0.3
v.number <- V(g1) |>length()
mean.length <- average.path.length(g1)
plot(g1, vertex.label = "", vertex.size = 2)
}
dev.off()


##-----------------------------------------------------------------------------
#Calculate path.length, k-edge and node degree for the possible targets.

list_path_conn_deg <- list()

for (i in seq_along(list_graphs)) {
  #All nodes
#gr.vert <-   t(combn(names(V(list_graphs[[i]])),2)) 
  #Only target nodes
gr.vert <- t(combn(possible_targets,2))
gr.vert_same <- cbind(possible_targets,possible_targets)
gr.vert <- rbind(gr.vert,gr.vert_same)

gr.vert <- gr.vert[,]
network <- list_graphs[[i]]
verti <- V(network) |> names()
ID <- 1:length(verti)
index <- cbind(verti,ID) |> as.data.frame() |> 
  mutate(ID=as.numeric(ID))

set <- gr.vert[,1]

set <- set[set %in% index$verti] |> as.data.frame()
colnames(set) <- "verti"

tab <- left_join(set,index,by="verti") |> 
  arrange(ID) |> distinct()
comb <- t(combn(tab$ID,2))


allcon <- future_apply(comb,1,function(edges){
  edge_connectivity(list_graphs[[i]],
                    source = edges[1], target = edges[2])
})
allpath <- future_apply(comb,1,function(edges){
  shortest.paths(list_graphs[[i]],
  v = edges[1], to = edges[2])
})
deg1 <- future_apply(comb,1,function(edges){
  igraph::degree(list_graphs[[i]],
                 v = edges[1])
})
deg2 <- future_apply(comb,1,function(edges){
  igraph::degree(list_graphs[[i]],
                 v = edges[2])
})



comb <- as.data.frame(comb)
names(comb) <- c("N1","N2")
comb$path.length <- allpath
comb$K.edge <- allcon
comb$Degree1 <- deg1
comb$Degree2 <- deg2

index1 <- index
colnames(index1) <- c("KEGG1","N1")
index2 <- index
colnames(index2) <- c("KEGG2","N2")

a <- left_join(comb,index1,by="N1")
b <- left_join(a,index2,by="N2")
comb <- b

comb <- comb |> 
  mutate(KEGG1_KEGG2 = paste0(KEGG1, "_", KEGG2)) |> 
  rowwise() |> 
  mutate(min_deg=min(Degree1,Degree2,na.rm=TRUE)) |> 
  mutate(max_deg=max(Degree1,Degree2,na.rm=TRUE)) |> 
  mutate(mean_deg=(Degree1+Degree2)/2) |> 
  ungroup() |> 
  distinct() |> 
  mutate(adjacency=ifelse(path.length==1,1,0))

list_path_conn_deg[[i]] <- comb
}

drug_mapping1 <- drug_mapping |>
rename(Drug1 = Drug,
KEGG1 = KEGG_eco)

df_join <- list_path_conn_deg |>
lapply(left_join, drug_mapping1, by = "KEGG1")

drug_mapping2 <- drug_mapping |>
rename(Drug2 = Drug,
KEGG2 = KEGG_eco)
df_join1 <- df_join |>
lapply(left_join, drug_mapping2, by = "KEGG2") |>
lapply(mutate, KEGG1_KEGG2 = paste0(KEGG1, "_", KEGG2)) |>
lapply(distinct)


df_join2 <- list()
for (j in seq_along(df_join1)) {
z <- df_join1[[j]]
for (i in 1:dim(z)[1]) {
a <- z[i, ]$Drug1
b <- z[i, ]$Drug2
ab <- sort(c(a, b))
z[i, ]$Drug1 <- ab[1]
z[i, ]$Drug2 <- ab[2]
}
df_join2[[j]] <- z
}

dist_conn_deg_adj <- df_join2 |>
lapply(mutate, drug_pair = paste0(Drug1, "_", Drug2)) |> 
lapply(select, KEGG1_KEGG2,
KEGG1, KEGG2, drug_pair,
Drug1, Drug2, path.length,
K.edge, Degree1, Degree2,adjacency,
mean_deg, min_deg, max_deg) |>
lapply(distinct) |>
lapply(filter, Drug1 != Drug2) |>
lapply(filter, Drug1 < Drug2) |> 
lapply(drop_na) |> 
lapply(arrange,drug_pair)

#Add rates and DDI info
rates <- read.csv(file.path(
  "data/4.PhylogeneticComparativeMethods/dataset1_ratios.csv"))

r <- rates |>
  select(clusters, sigma.rate) |>
  distinct() |>
  arrange(sigma.rate)
clusters <- factor(r$clusters, levels = r$clusters[order(r$sigma.rate)])

dist_conn_deg_adj <- lapply(dist_conn_deg_adj, inner_join, rates, by = "drug_pair") |>
  lapply(distinct)
names(dist_conn_deg_adj) <- net


df_target_tot  <- bind_rows(dist_conn_deg_adj, .id = "network") |>
  mutate(int_sign_ebw = factor(int_sign_ebw,
                               levels = c("Synergy", "Additivity", "Antagonism"))) |>
  mutate(drug_cat = ifelse(drug_category == "Same", 1, 0)) |>
  mutate(targ = ifelse(targeted_process == "Same", 1, 0)) |>
  mutate(us = ifelse(use == "Same", 1, 0)) |>
  mutate(ebw_g = ifelse(int_sign_ebw == "Antagonism",  +1,
                        ifelse(int_sign_ebw == "Additivity", 0,
                               ifelse(int_sign_ebw == "Synergy", -1, NA)))) |>
  mutate(ecr_g = ifelse(int_sign_ecr == "Antagonism",  +1,
                        ifelse(int_sign_ecr == "Additivity", 0,
                               ifelse(int_sign_ecr == "Synergy", -1, NA)))) |>
  mutate(seo_g = ifelse(int_sign_seo == "Antagonism",  +1,
                        ifelse(int_sign_seo == "Additivity", 0,
                               ifelse(int_sign_seo == "Synergy", -1, NA)))) |>
  mutate(stm_g = ifelse(int_sign_stm == "Antagonism",  +1,
                        ifelse(int_sign_stm == "Additivity", 0,
                               ifelse(int_sign_stm == "Synergy", -1, NA)))) |>
  mutate(pae_g = ifelse(int_sign_pae == "Antagonism",  +1,
                        ifelse(int_sign_pae == "Additivity", 0,
                               ifelse(int_sign_pae == "Synergy", -1, NA)))) |>
  mutate(pau_g = ifelse(int_sign_pau == "Antagonism",  +1,
                        ifelse(int_sign_pau == "Additivity", 0,
                               ifelse(int_sign_pau == "Synergy", -1, NA)))) |>
  rowwise() |>
  mutate(sum_g = ebw_g + ecr_g + seo_g + stm_g + pae_g + pau_g) |>
  distinct()|> 
  rename(Drug1=Drug1.x,Drug2=Drug2.x) |> 
  select(-c(Drug1.y,Drug2.y)) |> 
  distinct()

df_target_tot  |> 
  arrange(K.edge) |> 
  arrange(drug_pair) |> 
  arrange(network) |> 
write.csv(
  file.path("data/5.Targets_NetworkDistance",
            "df_target_tot_metrics.csv"), row.names = F)

significant_figures <- 4

final_DDI <- df_target_tot  |> 
  group_by(network,drug_pair) |> 
  summarise(mean.path.length = ifelse(mean(path.length,na.rm=TRUE)==Inf,Inf,round(mean(path.length,na.rm=TRUE),significant_figures)),
            mean.k.edge = round(mean(K.edge, na.rm = TRUE),significant_figures),
            mean.min.degree = round(mean(min_deg, na.rm = TRUE),significant_figures),
            mean.max.degree = round(mean(max_deg, na.rm = TRUE),significant_figures),
            mean.mean.degree = round(mean(mean_deg, na.rm = TRUE),significant_figures),
            max.adjacency = max(adjacency, na.rm = TRUE))

df1 <- df_target_tot |> select(drug_pair,Drug1,Drug2,code_3letter1,
                               code_3letter2,drug_category1,drug_category2,
                               targeted_cellular_process1,targeted_cellular_process2,
                               use1,use2,24:35,39:90)

df_DDI_tot  <- inner_join(final_DDI, df1, by = "drug_pair") |>
    distinct() |>
  mutate(max.adjacency = as.factor(max.adjacency)) |>
  mutate(connection_onoff =
           ifelse(mean.k.edge == 0, "disconnected", "connected")) 


df_DDI_tot   |> 
  arrange(mean.k.edge) |> 
  arrange(drug_pair) |> 
  arrange(network) |> 
  write.csv(
    file.path("data/5.Targets_NetworkDistance",
              "df_DDI_tot_metrics.csv"), row.names = F)

