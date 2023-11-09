# GO enrichment
library(pacman)
p_load(igraph,tidyverse,BioNet,ape,geomorph,matrixStats,
       ggpubr,plotrix,tidymodels,patchwork,here,EnvStats,report)

set_name <- "allddi" #allddi, sen2in1
# VISUALIZATION/STATS

df_target_tot <- read_csv(paste0("results/","df_target_metrics",".csv")) |> 
  mutate(KEGG1=toupper(KEGG1),
         KEGG2=toupper(KEGG2))

Uniprot <- read_csv(paste0("data/3.Targets_NetworkDistance/DrugTargets3_ecoli.csv")) 
  
Uniprot1 <- Uniprot |> 
  select(KEGG_eco,UNIPROT_eco) |> 
  distinct() |> 
  mutate(KEGG1=toupper(KEGG_eco)) |> 
  select(-KEGG_eco) |> 
  rename(Uniprot1=UNIPROT_eco)
Uniprot2 <- Uniprot |> 
  select(KEGG_eco,UNIPROT_eco) |> 
  distinct() |> 
  mutate(KEGG2=toupper(KEGG_eco))|> 
  select(-KEGG_eco) |> 
  rename(Uniprot2=UNIPROT_eco)

x1 <- left_join(df_target_tot,Uniprot1,by="KEGG1")
x2 <- left_join(x1,Uniprot2,by="KEGG2")


tab <-data.frame(cluster = numeric(0), 
                 targets = character(0))

all_clu <- x2 |>  
  filter(clusters_all%in% 1:16) |> 
  select(Uniprot1,Uniprot2) |> 
  pull() |> unique() |> sort()

t2 <- list()
for (i in sort(unique(x2$clusters_all))){
  t1 <- x2 |> 
    filter(clusters_all==i) |> 
    select(Uniprot1,Uniprot2) |> 
    pull() |> unique() |> sort()
  
  t2[[i]] <- paste0(t1,collapse = ",")
}

df <- do.call(rbind, t2) |> as.data.frame()
df$cluster <- 1:16


rownames(tab) <- NULL
colnames(tab) <- c("targets","clusters")

write.csv(tab,file = "results/allddi/targetsUNI_cluster.csv")

# Differential Enrichment

clu6 <- read_csv("results/allddi/enrichment/enrichment_clu6.csv") |> 
  filter(`Enrichment FDR`<=0.05) |> 
  mutate(GOID=str_extract(URL, "GO:\\d+")) |> 
  select(Pathway,GOID,nGenes)
clu13 <- read_csv("results/allddi/enrichment/enrichment_clu13.csv")|> 
  filter(`Enrichment FDR`<=0.05) |> 
  mutate(GOID=str_extract(URL, "GO:\\d+")) |> 
  select(Pathway,GOID,nGenes)
clu16 <- read_csv("results/allddi/enrichment/enrichment_clu16.csv")|> 
  filter(`Enrichment FDR`<=0.05) |> 
  mutate(GOID=str_extract(URL, "GO:\\d+")) |> 
  select(Pathway,GOID,nGenes)

library(clusterProfiler)
library(org.EcK12.eg.db)
library(UniprotR)

t2

split_list <- lapply(t2, function(x) unlist(strsplit(x, ",")))

result <- lapply(t2, function(vec) ConvertID(vec, 
                                             ID_from = "UniProtKB_AC-ID", 
                                             ID_to = "GeneID"))
result <- lapply(result, unname)

gene_list_split <- map(result, ~ str_split(.x, ",")) |> 
  unlist(recursive = FALSE)

ck <- compareCluster(geneCluster = gene_list_split, 
                     fun = enrichKEGG)

ggo <- groupGO(gene     = gene,
               OrgDb    = org.EcK12.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)

