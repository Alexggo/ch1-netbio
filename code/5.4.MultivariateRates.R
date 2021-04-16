library(geomorph)
library(tidyverse)
library(broom)
library(OUwie)

mat <- read.csv("data/processed/Broc2018_maindataset.csv")


for (i in 1:dim(mat)[1]){
  drug1 <- mat$drug1[i]
  drug2 <- mat$drug2[i]
  vector1 <- sort(c(drug1,drug2))
  mat$drugID[i] <- paste0(vector1[1],"_",vector1[2])
}

mat$test <- mat$drug_pair==mat$drugID
mat %>% select(drug_pair,drugID,test) %>% View()




  filter(!is.na(ebw)&!is.na(ecr)&!is.na(seo)&!is.na(stm)&!is.na(pae)&!is.na(pau))




clust <- read.csv("results/tSNE_Clustermembership_p706.csv")
colnames(clust)[1] <- "drug_pair"

mat1 <- full_join(mat,clust,by="drug_pair")

# Read phylogenetic tree
species <- c("ebw","ecr","stm","seo","pae","pau")

treetype <- "bac_species.txt.hcp"
treefile <- file.path('data/PhySpeTree',treetype,"iqtree.tree.treefile")
tree_nw <- read.tree(treefile)
# Drop the species that are not needed.
not.species <- tree_nw$tip.label[!(tree_nw$tip.label %in% species)]
tree_nw <- drop.tip(tree_nw,not.species)

a <- mat1 %>% as_tibble() %>% group_by(clusters) %>%
  summarize(mean.ebw=mean(ebw),
            mean.ecr=mean(ecr),
            mean.seo=mean(seo),
            mean.stm=mean(stm),
            mean.pae=mean(pae),
            mean.pau=mean(pau))

# Each cluster of traits is a "trait".

b <- a %>% select(-1) %>% t()
rownames(b) <- c("ebw","ecr","seo","stm","pae","pau")
colnames(b) <- paste0("Cluster-",LETTERS[1:13])



EvoModels <- c("BM1","OU1")
fitEvoModels <- list()
AICc <- list()
BM1table <- list()
OU1table <- list()
d1_full <- b

d <- d1_full+1 # I add 1 to all values so that they are not negative.

d <- d #Some DDI give an error and have to be excluded or the loop stops.


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
  OU1table[[i]] <- fitEvoModels[[i]][[2]]$solution
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
  mutate(AICc.OU.BM=OU1-BM1) %>%
  arrange(AICc.OU.BM,decreasing=T)


ggplot(df.range) +
  aes(x = drug_pair, fill = min.model, weight = AICc.OU.BM) +
  geom_bar() +
  scale_fill_hue() +
  theme_minimal()

