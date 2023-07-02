#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#' ---
#' title: "4.Heatmaps_PCA_BlombergK"
#' author: "Alex Gil"
#' date: "1/26/2021"
#' output: html_document
#' 
#' # 4.ChemicalMatrices.
#' 
#' 
## -----------------------------------------------------------------------------
library("ape")
library("phangorn")
library('phytools')
library("treeio")
library("ggtree")
library('tidyverse')
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
#' 
## -----------------------------------------------------------------------------

# Read data matrix
matrixfile <- file.path("data/processed","2018_Broc_ED02.1.csv")

mat_traits <- read.csv(matrixfile,row.names = 1) %>% t()

#Remove columns with at least one NA.
mat_traits <- mat_traits[ , apply(mat_traits, 2, function(x) !any(is.na(x)))]
#' ## Remove zero variance columns from mat_traits.
## -----------------------------------------------------------------------------
d1_full <- mat_traits[ , colSums(is.na(mat_traits)) == 0]
which(apply(d1_full, 2, var)==0)
d1_full <- d1_full[ , which(apply(d1_full, 2, var) != 0)]

species <- mat_traits %>% rownames()

# Read phylogenetic tree
treetype <- "bac_species.txt.hcp"
treefile <- file.path('data/PhySpeTree',treetype,"iqtree.tree.treefile")
tree_nw <- read.tree(treefile)

# Drop the species that are not needed.
not.species <- tree_nw$tip.label[!(tree_nw$tip.label %in% species)]
tree_nw <- drop.tip(tree_nw,not.species)

#pdf(file.path("results",paste0(args[2],".pdf")))

#Plot tree
tree <- ggtree(tree_nw,layout = 'rectangular')+
  geom_tiplab(hjust=-0.5,fontface="bold",align = T, color='steelblue')+
  geom_treescale(fontsize = 3)

tree

#' ## Calculate dendogram for traits (heatmap) and plot tree.
## -----------------------------------------------------------------------------
dend <- hclust(dist(mat_traits,method="euclidean"),method="average")  

dend %>% 
  as.phylo() %>%   ggtree(layout = 'rectangular')+
  geom_tiplab(hjust=-0.5,fontface="bold",align = T, color='steelblue')+
  geom_treescale(fontsize = 3)

Heatmap(mat_traits,show_column_names = F,
        row_names_side = "left",
        column_dend_height = unit(5,"cm"),
        cluster_rows = color_branches(dend,k=3),
        show_column_dend = T,row_title = "Strains",
        column_title_side = "bottom",column_title = "Pairwise drug-drug interaction",
        row_names_gp = gpar(cex=0.75),row_dend_width = unit(2,"cm"),
        name = "Score",col=rev(brewer.pal(10,'RdBu')))


## Plot PCAs, and phylomorphospaces using d1_full matrix.
## -----------------------------------------------------------------------------
pca=prcomp(d1_full, center = TRUE, scale. = TRUE)
pca$names=row.names(d1_full)
## Eigenvalues
fviz_eig(pca)
## PCA plot
fviz_pca_ind(pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)     # Avoid text overlapping

#PCA
summary(pca)
pca.x=pca$x
pca.1.2=pca.x[,c(1,2)]

#Phylomorphospace plots.
phylomorphospace(tree_nw, pca$x[,1:2], label = c("horizontal"), node.size=c(.5,1))
phylomorphospace3d(tree_nw,pca.x[,c(1,2,3)],method="static")
fancyTree(tree_nw,X=pca.x[,c(1,2)],type="traitgram3d",method="static")

tree_log <- tree_nw
tree_log$edge.length <- log(tree_nw$edge.length)-min(log(tree_nw$edge.length))+1
fancyTree(tree_log,X=pca.x[,c(1,2)],type="traitgram3d",method="static")


## Eigenvectors plots.
fviz_pca_biplot(pca, 
                repel = TRUE,
                col.var = "contrib", # Color by contributions to the PC
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                col.ind = "#696969")  # Individuals color




# Which drug combinations have the highest sd?
x <- apply(d1_full,2,sd) %>%
  sort(decreasing = T)

par(mfrow=c(1,1))
highest_var<- x[1:4]
mred <- d1_full[,names(highest_var)]
obj<-fancyTree(tree_nw,type="scattergram",X=mred)
plot(obj)

par(mfrow=c(2,2))

for (i in 1:dim(mred)[2]){
  fancyTree(tree_nw,type="phenogram95",x=mred[,i],spread.cost=c(1,0))
  title(colnames(mred)[i])
}

for (i in 1:dim(mred)[2]){
  obj<-contMap(tree_nw,mred[,i],plot=FALSE)
  obj<-setMap(obj,invert=TRUE)
  n <- colnames(mred)[i]
  plot(obj,fsize=c(0.4,1),outline=FALSE,lwd=c(3,7),main=as.character(n))
  title(colnames(mred)[i])
}

m1 <- mat_traits[,names(highest_var)]
# Dendrograms from matrices
trees <- list()
for (i in 1:4){
  trees[[i]] <- hclust(dist(m1[,i],method="euclidean"),method="average")  %>% 
    as.phylo
}

par(mfrow=c(2,2))

for (i in 1:4){
  trees[[i]] %>% plot()
}

for (i in 1:dim(mred)[2]){
  obj<-contMap(trees[[i]],mred[,i],plot=FALSE)
  obj<-setMap(obj,invert=TRUE)
  n <- colnames(mred)[i]
  plot(obj,fsize=c(0.4,1),outline=FALSE,lwd=c(3,7),main=as.character(n))
  title(colnames(mred)[i])
}

par(mfrow=c(1,1))

#' Which are the DDIs with highest Phylogenetic signal?
## -----------------------------------------------------------------------------
# Blomberg’s K for each DDI.
BlombergK <- apply(d1_full,2,phylosig,tree=tree_nw,method = "K", test = TRUE)


K <-sapply(BlombergK, function(x)x[1]) %>% unlist()
P <-sapply(BlombergK, function(x)x[2]) %>% unlist()

df.K <- data.frame(Kappa=K,
                   P.value=P,
                   DDI=names(K))  %>% 
  as_tibble() %>%
  mutate(DDI=str_sub(DDI, end=-3)) %>% 
  mutate(phylo.DDI=ifelse(P.value>0.05,"Non-significant",ifelse(Kappa>0.3,"High and Significant","Low and Significant"))) %>% 
  mutate(If.DDI=ifelse(phylo.DDI=="High and Significant",DDI,""),
         If.test=ifelse(phylo.DDI=="High and Significant",T,F))

df.K %>% 
  ggplot(aes(x=P.value,y=Kappa,col=phylo.DDI,label=If.DDI))+
  theme_minimal()+
  geom_vline(xintercept = 0.05)+
  labs(x="P value",
       y="Phylogenetic Signal (Blomberg's Kappa)",
       title="Phylogenetic Signal for each DDI",
       subtitle=paste("File:",matrixfile),
       caption = "K close to zero correspond to a random or convergent pattern of evolution.
       K = 1, corresponds to BM. 
       K>1 indicates strong phylogenetic signal.")+
  geom_point()

 

df.K %>% 
  filter(P.value<=0.06) %>% 
  ggplot(aes(x=P.value,y=Kappa,col=phylo.DDI,label=If.DDI))+
  theme_minimal()+
  geom_vline(xintercept = 0.05,linetype="dashed")+
  labs(x="P value",
       y="Phylogenetic Signal (Blomberg's Kappa)",
       title="Phylogenetic Signal for each DDI",
       subtitle=paste("File:",matrixfile),
       caption = "K close to zero correspond to a random or convergent pattern of evolution.
       K = 1, corresponds to BM. 
       K>1 indicates strong phylogenetic signal.")+
  geom_point()+
  geom_text(aes(label=ifelse(P.value<0.05 & Kappa >0.35,as.character(DDI),'')),hjust=-0.05,vjust=0)





