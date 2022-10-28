library(pacman)
p_load(bigMap,bigmemory,tidyverse,rlist)


# Dataset
file <- file.path("data/1.processed/Broc2018_maindataset.csv")
dataset <- read_csv(file) %>%
  filter(!is.na(ebw)&!is.na(ecr)&!is.na(seo)&!is.na(stm)&!is.na(pae)&!is.na(pau))


# Explore a list of perplexities

files <- list.files(path = "results",pattern = "\\.RData$")
files <- files[c(9,7,10,11,1:6)]
files <- paste0("results/",files)

# Load lists
load(files[1])
list1 <- m.list1
load(files[2])
list2 <- m.list1
load(files[3])
list3 <- m.list1
load(files[4])
list4 <- m.list1
load(files[5])
list5 <- m.list1
load(files[6])
list6 <- m.list1
load(files[7])
list7 <- m.list1
load(files[8])
list8 <- m.list1
load(files[9])
list9 <- m.list1
load(files[10])
list10 <- m.list1

big_list <- c(list1,list2,list3,list4,list5,
              list6,list7,list8,list9,list10)
ptsne.ppx <- seq(5,2505,10)

pdf()
for (i in 1:length(big_list)){
  bdm.cost(big_list[[i]])
  title(main = paste("perplexity = ",as.character(ptsne.ppx[i])))
  bdm.ptsne.plot(big_list[[i]], ptsne.cex = 4)
  title(main = paste("perplexity = ",as.character(ptsne.ppx[i])))  
  bdm.wtt.plot(big_list[[i]])
  title(main = paste("perplexity = ",as.character(ptsne.ppx[i])))
}
dev.off()

x <- 1:251
index <- x[ptsne.ppx==805]
# Find a value that makes the size function stable.
bdm.cost(big_list[[index]])
bdm.ptsne.plot(big_list[[index]], ptsne.cex = 2)
bdm.wtt.plot(big_list[[index]])
# The embedding function appears to reach stability the higher the ppx.
# However, for higher ppx, fewer clusters are formed, but the variability between threads increases.
# Cost and effect size curves look reasonable for ppx above 235. 2255 if the largest ppx with 2 clusters.


# Get cluster labels
list_clusters <- list()
for (i in 1:length(big_list)){
  list_clusters[[i]] <- bdm.labels(big_list[[i]])
}

x <- do.call(rbind, list_clusters)
rownames(x) <- paste0("ppx_",ptsne.ppx)
y <- t(x) |> as.data.frame()
y$drugdrug <- dataset$drugdrug

full_dataset <- left_join(dataset,y,by="drugdrug") |> 
  filter(!is.na(ebw)&!is.na(ecr)&!is.na(seo)&!is.na(stm)&!is.na(pae)&!is.na(pau))


## Plot labels in the graph.
df <- dataset%>%mutate(across(where(is.character), as.factor)) |> 
  mutate(across(where(is.factor), as.numeric))




category <- df |> select(1:86) |> as.matrix()
bdm.qMap(big_list[[1]], data = category, qMap.cex = 1)

write.csv(full_dataset,"data/3.InteractionScores_tSNE/all_ppx_large.csv",row.names = F)
