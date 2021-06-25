#' ---
#' title: "4.2.OUwie_pertrait.R"
#' author: "Alex Gil"
#' date: "1/26/2021"
#' output: html_document
#' ---
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
library("ComplexHeatmap")
library("dendextend")
library("RColorBrewer")
library("car")
library("factoextra")
library("ggbiplot")
library("OUwie")
library("ouch")
library("factoextra")
library("geiger")
library("ggrepel")
library("hrbrthemes")
library("ggforce")
library("patchwork")

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
#'     + 2018_Broc_ED09C.csv
#' 
## -----------------------------------------------------------------------------

# Read data matrix
file <- file.path("data/1.processed","Broc2018_maindataset.csv")

#Remove DDI with at least one NA.
dataset <- read_csv(file) %>%
  filter(!is.na(ebw)&!is.na(ecr)&!is.na(seo)&!is.na(stm)&!is.na(pae)&!is.na(pau))

mat_red <- dataset %>% select(ebw,ecr,stm,seo,pae,pau) %>%as.matrix
rownames(mat_red) <- dataset$drug_pair

mat_red <- mat_red %>% t()

#Remove columns with at least one NA.
mat_red <- mat_red[ , apply(mat_red, 2, function(x) !any(is.na(x)))]
#' ## Remove zero variance columns from mat_traits.
## -----------------------------------------------------------------------------
mat_red <- mat_red[ , colSums(is.na(mat_red)) == 0]
mat_red <- mat_red[ , which(apply(mat_red, 2, var) != 0)]

species <- mat_traits %>% rownames()

# Read phylogenetic tree
treetype <- "bac_species.txt.hcp"
treefile <- file.path('data/2.Phylogenetics_PhySpeTree',treetype,"iqtree.tree.treefile")
tree_nw <- read.tree(treefile)

# Drop the species that are not needed.
not.species <- tree_nw$tip.label[!(tree_nw$tip.label %in% species)]
tree_nw <- drop.tip(tree_nw,not.species)

#Tree should have internal node labels (for regimes).
tree_nw$node.label <- as.character(c(1,1,1,1,1))

#Plot tree
tree <- ggtree(tree_nw,layout = 'rectangular')+
  geom_tiplab(hjust=-0.5,fontface="bold",align = T, color='steelblue')+
  geom_treescale(fontsize = 3)

tree

## Model traits using OUwie.
# BM1
# BMS
# OU1
# OUM
# alpha and sigma^2 vary (i.e. OUMV, OUMA, and OUMVA)
# OUMV
# OUMV
# OUMVA

EvoModels <- c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA")
fitEvoModels <- list()
AICc <- list()
BM1table <- list()
OU1table <- list()
d1_full <- mat_red

d <- d1_full+1 # I add 1 to all values so that they are not negative.

d <- d[,-c(128,406,486,1282,1365,1518,2234,2592)] #Some DDI give an error and have to be excluded or the loop stops.


for (i in 1:dim(d)[2]){
  trait <- data.frame(Genus_species=rownames(d),
             Reg=1,
             X=d[,i])
  print(colnames(d)[i])
  print(i)
  print(d[,i])
  fitEvoModels[[i]] <- lapply(EvoModels,
                              FUN = function(x) OUwie(tree_nw, trait, model=x, root.station=FALSE, 
                                                      scaleHeight=TRUE, algorithm="invert", quiet=TRUE,
                                                      simmap.tree = FALSE,check.identify = FALSE))
  AICc[[i]] <- sapply(fitEvoModels[[i]], 
                      FUN = function(x) x$AICc)
  names(AICc[[i]]) <- EvoModels
  BM1table[[i]] <- fitEvoModels[[i]][[1]]$solution
  OU1table[[i]] <- fitEvoModels[[i]][[3]]$solution
}

# For AICc
table <- as.data.frame(do.call(rbind, AICc))
vect <- d %>% colnames()
table$drug_pair <- vect

table1 <- table %>% as.tibble() %>% 
  mutate(min.model=colnames(table)[apply(table,1,which.min)])

# Parameters for BM
df.BM <- as.data.frame(do.call(rbind, BM1table))
df.BM$parameter <- c("alpha","sigma.sq")
df.BM$drug_pair <- rep(table$drug_pair,each=2)
colnames(df.BM) <- c("value","parameter","drug_pair")
rownames(df.BM) <- NULL
df.BM <- df.BM %>% select(drug_pair,parameter,value) %>% 
  pivot_wider(values_from = value,names_from=parameter)%>% 
  mutate(BM.alpha=alpha,BM.sigmasq=sigma.sq) %>%
  select(drug_pair,BM.alpha,BM.sigmasq)
# Parameters for OU
df.OU <- as.data.frame(do.call(rbind, OU1table))
df.OU$parameter <- c("alpha","sigma.sq")
df.OU$drug_pair <- rep(table$drug_pair,each=2)
colnames(df.OU) <- c("value","parameter","drug_pair")
rownames(df.OU) <- NULL
df.OU <- df.OU %>% select(drug_pair,parameter,value)%>% 
  pivot_wider(values_from = value,names_from=parameter) %>% 
  mutate(OU.alpha=alpha,OU.sigmasq=sigma.sq) %>%
  select(drug_pair,OU.alpha,OU.sigmasq)

df.all <- full_join(df.BM,df.OU)
df.all <- full_join(table1,df.all) %>%
  select(drug_pair,BM1,OU1,min.model,BM.sigmasq,OU.alpha,OU.sigmasq)
df.range <- df.all %>%
  mutate(AICc.range=abs(OU1-BM1)) %>%
  arrange(AICc.range,decreasing=T)
write.csv(df.range,"data/5.PhylogeneticComparativeMethods/ED09C_tableAICc.csv",row.names = F)






df.red <- df.range %>% filter(drug_pair %in% colnames(mat_red))

# Get the 4 highest-AICc difference OU1, and the 4 highest-AICc difference BM1.
OU <- df.red %>% filter(min.model=="OU1") %>% 
  arrange(AICc.range,decreasing=T) %>% filter(AICc.range>33.5)
BM <- df.red %>% filter(min.model=="BM1") %>% 
  arrange(AICc.range,decreasing=T) %>% filter(AICc.range>8.7)

test0 <- rbind(OU,BM) %>% as.data.frame() 
test <- test0[,1]

# First four are OU1, last 4 are BM1
m1 <- mat_traits[,test]

# Plot trees of the DDI with highest difference in AICcs between models.


#------------------
p <- list()
for (i in 1:8){
  svl <- as.matrix(m1)[,i]
  fit <- phytools::fastAnc(tree_nw,svl,vars=TRUE,CI=TRUE)
  
  td <- data.frame(node = nodeid(tree_nw, names(svl)),
                   trait = svl)
  nd <- data.frame(node = names(fit$ace), trait = fit$ace)
  
  d <- rbind(td, nd)
  d$node <- as.numeric(d$node)
  tree <- full_join(tree_nw, d, by = 'node')
  
  p[[i]] <- ggtree(tree, aes(color=trait), continuous = TRUE, yscale = "trait") + 
    scale_color_viridis_c() + theme_minimal() +geom_tiplab(col="black",size=3)+
    ggtitle(colnames(m1)[i],subtitle = paste0(test0$min.model[i],". OU Alpha=",round(test0$OU.alpha[i],2))) #+ ylim(-0.4,0.1)
}  


(p[[1]]+p[[2]]+p[[3]]+p[[4]])/(p[[5]]+p[[6]]+p[[7]]+p[[8]])

trees <- list()
for (i in 1:8){
  trees[[i]] <- hclust(dist(m1[,i],method="euclidean"),method="average")  %>% 
    as.phylo
}

par(mfrow=c(2,4))

for (i in 1:8){
  trees[[i]] %>% plot()
}

for (i in 1:8){
  obj<-contMap(tree_nw,m1[,i],plot=FALSE)
  obj<-setMap(obj,invert=TRUE)
  n <- colnames(m1)[i]
  plot(obj,fsize=c(0.4,1),outline=FALSE,lwd=c(3,7),main=as.character(n))
  title(colnames(m1)[i])
}

for (i in 1:8){
  fancyTree(tree_nw,type="phenogram95",x=m1[,i],spread.cost=c(1,0))
  title(colnames(m1)[i])
}
