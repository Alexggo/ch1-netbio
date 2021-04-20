#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#' ---
#' title: 3.3.DistanceFromMatrix.
#' author: "Alex Gil"
#' 
#' 
## -----------------------------------------------------------------------------
library("plyr")
library("tidyverse")
library("broom")
library("ape")
library("phangorn")
library('phytools')
library("treeio")
library("ggtree")
library('plotly')
library("dendextend")
library("RColorBrewer")
library("car")
library("factoextra")
library("OUwie")
library("ouch")
library("factoextra")
library("geiger")
library("ggrepel")
library("hrbrthemes")
library("ggforce")
library("ComplexHeatmap")
library("ggbiplot")


#' Input files, phylogenetic tree and matrix.
#' 
#'   -Choose the tree: args[1]
#'     + bac_species.txt.hcp
#'     + bac_species.txt.rna
#'     + fun_species.txt.hcp
#'     + fun_species.txt.rna
#'   -Choose the matrix (rows are species, columns are traits): args[2]
#'     + 2011_Spit_S1.1.csv
#'     + 2011_Spit_S3.1.csv
#'     + 2015_Spit_S2.1.csv
#'     + 2015_Spit_S5.1.csv
#'     + 2018_Broc_ED02.1.csv
#'     + 2018_Broc_ED03.2.csv
#'     + 2018_Broc_ED09C.csv.mat.csv
#' 
## -----------------------------------------------------------------------------

mat <- read.csv("data/1.processed/Broc2018_maindataset.csv") %>%
  filter(!is.na(ebw)&!is.na(ecr)&!is.na(seo)&!is.na(stm)&!is.na(pae)&!is.na(pau))

clust <- read.csv("data/4.InteractionScores_tSNE/tSNE_Clustermembership_ppx706.csv")

mat1 <- full_join(mat,clust,by="drugdrug")

# Read phylogenetic tree
species <- c("ebw","ecr","stm","seo","pae","pau")

treetype <- "bac_species.txt.hcp"
treefile <- file.path('data/2.Phylogenetics_PhySpeTree',treetype,"iqtree.tree.treefile")
tree_nw <- read.tree(treefile)
# Drop the species that are not needed.
not.species <- tree_nw$tip.label[!(tree_nw$tip.label %in% species)]
tree_nw <- drop.tip(tree_nw,not.species)

phylofile <- file.path("data/2.Phylogenetics_PhySpeTree","dist_bac_hcp.csv")
PD_dist <- read.csv(file=phylofile)%>%
  mutate(ID=paste0(Sp1,"-",Sp2)) %>% 
  filter(PhyloDist!=0) 


# Dendrograms from matrices

grouped <- mat1 %>% group_by(clusters) %>% select(12:17)

list1 <- list()
for (i in 1:length(unique(grouped$clusters))){
  m <- grouped %>% filter(clusters==i) %>% ungroup() %>% 
    select(-clusters) %>% t()
  
  dend <- hclust(dist(m,method="euclidean"),method="average")  %>% 
    as.phylo()
  
  PatristicDistMatrix<-cophenetic(dend)
  tabnew<- PatristicDistMatrix %>% as.data.frame()
  
  tabnew <- tabnew[ , order(names(tabnew))]
  tabnew$Sp1 <- rownames(tabnew)
  
  tabnew <- tabnew %>%
    arrange(Sp1) %>% 
    select(-Sp1)
  
  tabnew[lower.tri(tabnew)] <- NA
  tabnew$Sp1 <- colnames(tabnew)
  
  tabnew <- tabnew %>%
    select(dim(tabnew)[2],1:(dim(tabnew)[2]-1))
  
  tabnew <- tabnew  %>% 
    pivot_longer(names_to="Sp2",values_to="Chemdist",2:(dim(tabnew)[2]))%>%
    filter(!is.na(Chemdist)) %>% 
    filter(Chemdist!=0) %>% 
    mutate(cluster=paste0("Cluster_",i),
           ID=paste0(Sp1,"-",Sp2))
  
  list1[[i]] <- tabnew
}
Chem_Dist <- as.data.frame(do.call(rbind, list1))

#Merge datasets.
df <- inner_join(PD_dist,Chem_Dist,by="ID")  %>% 
  select(!c(Sp1.x,Sp2.x)) %>% 
  rename(Sp1=Sp1.y,Sp2=Sp2.y) %>% 
  filter(!is.na(Chemdist)) %>% 
  mutate(S1 = ifelse(Sp1 %in% c("ebw","ecr"),"E. coli",ifelse(Sp1 %in% c("pae","pau"),"P. aeruginosa",ifelse(Sp1 %in% c("seo","stm"),"S. enterica",ifelse(Sp1 =="cal","C. albicans",ifelse(Sp1 =="cne","C. neoformans",ifelse(Sp1=="sce","S. cerevisiae",ifelse(Sp1=="cgi","C. gattii",ifelse(Sp1=="spo","S. pombe","Unknown")))))))))%>% 
  mutate(S2 = ifelse(Sp2 %in% c("ebw","ecr"),"E. coli",ifelse(Sp2 %in% c("pae","pau"),"P. aeruginosa",ifelse(Sp2 %in% c("seo","stm"),"S. enterica",ifelse(Sp2 =="cal","C. albicans",ifelse(Sp2 =="cne","C. neoformans",ifelse(Sp2=="sce","S. cerevisiae",ifelse(Sp2=="cgi","C. gattii",ifelse(Sp2=="spo","S. pombe","Unknown"))))))))) %>% 
  mutate(class=ifelse(S1==S2,"Intraspecific comparison",ifelse(S1=="P. aeruginosa","Interspecific comparison (with Pseudomonas)",ifelse(S2=="P. aeruginosa","Interspecific comparison (with Pseudomonas)","Interspecific comparison (other than Pseudomonas)")))) %>% 
  mutate(cluster=fct_relevel(cluster,levels=paste0("Cluster_",1:13)))


#Model statistics.
lmod <- lm(Chemdist~log(PhyloDist),data=df)
tid <- tidy(lmod)
aug <- augment(lmod)
gla <- glance(lmod)

ggplot(df,aes(x=PhyloDist,y=Chemdist,label=ID,group=cluster))+
  geom_point(aes(color=class))+
  geom_text_repel()+
  geom_smooth(method = "lm", formula=y~log(x),se=FALSE, color="black")+
  labs(x="Phylogenetic Distance",
       y="Chemical Patristic Distance from Dendrogram",
       title="Chemical divergence accumulates with time")+
  theme_ipsum_rc(grid="XY")+
  theme_minimal()+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  facet_wrap(.~cluster)



ggplot(df,aes(x=PhyloDist,y=Chemdist,label=ID))+
  geom_point(aes(color=class))+
  geom_text_repel()+
  geom_smooth(method = "lm", se=FALSE, color="black")+
  labs(x="Phylogenetic Distance",
       y="Chemical Patristic Distance from Dendrogram",
       title="Chemical divergence accumulates with time")+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  annotation_logticks(sides="b")+
  theme_ipsum_rc(grid="XY")+
  theme_minimal()+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  facet_wrap(.~cluster)

