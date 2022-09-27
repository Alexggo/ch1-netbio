library(pacman)
p_load(bigMap,bigmemory,tidyverse,rlist)

# Dataset
file <- file.path("data/1.processed/Broc2018_maindataset.csv")
dataset <- read_csv(file) %>%
  filter(!is.na(ebw)&!is.na(ecr)&!is.na(seo)&!is.na(stm)&!is.na(pae)&!is.na(pau)) %>% 
  mutate(useuse=as.factor(useuse),
         categorycategory=as.factor(categorycategory),
         processprocess=as.factor(processprocess),
         drug_category=as.factor(drug_category),
         use=as.factor(use),
         targeted_process=as.factor(targeted_process),
         int_sign_ebw=as.factor(int_sign_ebw))


# Explore a large range of perplexities.
load('results/B.SameRes_across_strains_DDIs/RObjects/mlist1.RData')
list1 <- m.list1
ptsne.ppx <- c(seq(5,50,10),seq(50, 2120, by = 82))

pdf("img/B.SameRes_across_strains_DDIs/tSNE.pdf")
for (i in 1:length(list1)){
  bdm.cost(list1[[i]])
  title(main = paste("perplexity = ",as.character(ptsne.ppx[i])))
  bdm.ptsne.plot(list1[[i]], ptsne.cex = 4)
  title(main = paste("perplexity = ",as.character(ptsne.ppx[i])))  
  bdm.wtt.plot(list1[[i]])
  title(main = paste("perplexity = ",as.character(ptsne.ppx[i])))
}
dev.off()

# Perplexity ranges between 5 and 2100.
#5   15   25   35   45   50  132  214  296  378  460  542  624  706  788  870  952 1034 1116
# 1198 1280 1362 1444 1526 1608 1690 1772 1854 1936 2018 2100
x <- 1:31
index <- x[ptsne.ppx==132]
# Find a value that makes the size function stable.
bdm.cost(list1[[index]])
bdm.ptsne.plot(list1[[index]], ptsne.cex = 2)
bdm.wtt.plot(list1[[index]])
# With a ppx above 706 the embedding function appears to reach stability.
# For ppx>=870, only 5 clusters are formed, but the variability between threads increases.


## Plot labels in the graph.
labels <- dataset %>% 
  select(drug_pair) %>% pull() # Select your clustering label processprocess, useuse, categorycategory.
ID <- dataset %>% 
  select(drugdrug) %>% pull()

subset <- as.factor(labels) %>% unique()



### Looking at the cluster membership
Cluster_membership <- data.frame(drugdrug=ID,
                                 drug_pair=labels, 
                                 clusters=bdm.labels(list1[[index]]))



dat <- as.matrix(dataset[,12:17])
bdm.qMap(list1[[index]], data = dat, qMap.cex = 1)

# Add labels to plot
classG<- as.numeric(labels)
list1[[index]]$lbls <- classG


pdf()
for (i in 1:length(subset)){
  bdm.qMap(list1[[index]], 
           data = list1[[index]]$data[, 1:6], 
           subset = which(list1[[index]]$lbls %in% i))
  title(sub = levels(subset)[i])
}
dev.off()


write.csv(Cluster_membership,"data/3.InteractionScores_tSNE/tSNE_Clustermembership_ppx706.csv",row.names = F)



# Explore a smaller range of perplexities closer to optimal value
load('results/RObjects/mlist2_ppx0700_0900.RData')
list2 <- m.list1
ptsne.ppx <- ptsne.ppx <- seq(50, 600, by = 10)

pdf()
for (i in 1:length(list2)){
  bdm.cost(list2[[i]])
  title(main = paste("perplexity = ",as.character(ptsne.ppx[i])))
  bdm.ptsne.plot(list2[[i]], ptsne.cex = 4)
  title(main = paste("perplexity = ",as.character(ptsne.ppx[i])))  
  bdm.wtt.plot(list2[[i]])
  title(main = paste("perplexity = ",as.character(ptsne.ppx[i])))
}
dev.off()


# Perplexity ranges between 5 and 2100.
#700 725 750 775 800 825 850 875 900
x <- 1:9
index <- x[ptsne.ppx==825]
# Find a value that makes the size function stable.
bdm.cost(list2[[index]])
bdm.ptsne.plot(list2[[index]], ptsne.cex = 2)
bdm.wtt.plot(list2[[index]])
# For ppx=825, cost function reaches stable solution.



## Plot labels in the graph.
labels <- dataset %>% 
  select(categorycategory) %>% pull() # Select your clustering label processprocess, useuse, categorycategory.
ID <- dataset %>% 
  select(drugdrug) %>% pull()

subset <- as.factor(labels) %>% unique()




### Looking at the cluster membership
Cluster_membership <- data.frame(drugdrug=ID,
                                 drug_pair=labels, 
                                 clusters=bdm.labels(list2[[index]]))



dat <- as.matrix(dataset[,12:17])
bdm.qMap(list2[[index]], data = dat, qMap.cex = 1)

# Add labels to plot
classG<- as.numeric(labels)
list2[[index]]$lbls <- classG


pdf()
for (i in 1:length(subset)){
  bdm.qMap(list2[[index]], 
           data = list2[[index]]$data[, 1:6], 
           subset = which(list2[[index]]$lbls %in% i))
  title(sub = levels(subset)[i])
}
dev.off()

write.csv(Cluster_membership,"data/3.InteractionScores_tSNE/tSNE_Clustermembership_ppx825.csv",row.names = F)

