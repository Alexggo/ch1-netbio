library(pacman)
p_load(bigMap,bigmemory,tidyverse,rlist)

set_name <- "onlysen" #allddi,sen2in1,onlysen
ptsne.ppx <- seq(5,2505,10)

# Merge tsne
# load(file="results/set_name/tSNE_5_255.RData")
# f1 <- m.list1
# load(file="results/set_name/tSNE_265_505.RData")
# f2 <- m.list1
# load(file="results/set_name/tSNE_515_755.RData")
# f3 <- m.list1
# load(file="results/set_name/tSNE_765_1005.RData")
# f4 <- m.list1
# load(file="results/set_name/tSNE_1015_1255.RData")
# f5 <- m.list1
# load(file="results/set_name/tSNE_1265_1505.RData")
# f6 <- m.list1
# load(file="results/set_name/tSNE_1515_1755.RData")
# f7 <- m.list1
# load(file="results/set_name/tSNE_1765_2005.RData")
# f8 <- m.list1
# load(file="results/set_name/tSNE_2015_2255.RData")
# f9 <- m.list1
# load(file="results/set_name/tSNE_2265_2505.RData")
# f10 <- m.list1
# big_list <- c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10)
# save(big_list,file="results/set_name/all_tsne.RData")


# Explore a list of perplexities
# Load lists
load(file=paste0("results/",set_name,"/all_tsne.RData"))
ptsne.ppx <- seq(5,2505,10)

#Make figure

pdf(paste0("results/",set_name,"/",set_name,".pdf"))

for (i in 1:length(big_list)){
  bdm.cost(big_list[[i]])
  title(main = paste("perplexity = ",as.character(ptsne.ppx[i])))

  bdm.ptsne.plot(big_list[[i]], ptsne.cex = 4)
  title(main = paste("perplexity = ",as.character(ptsne.ppx[i])))

  
  bdm.wtt.plot(big_list[[i]])
  title(main = paste("perplexity = ",as.character(ptsne.ppx[i])))

}

dev.off()


ptsne.ppx <- seq(5,2505,10)
# ALL DDI. Good tsnes between 205-1455. Indexes 21-146
load(file=paste0("results/allddi/all_tsne.RData"))
allddi <- big_list
# SET2. Good tsnes between 115-455. Indexes 12-46
load(file=paste0("results/sen2in1/all_tsne.RData"))
sen2in1 <- big_list
# ONLYSEN. Good tsnes between 145-455. Index 15-46
load(file=paste0("results/onlysen/all_tsne.RData"))
onlysen <- big_list

list_ppx <- list(allddi,sen2in1,onlysen)
names(list_ppx) <- c("allddi","sen2in1","onlysen")
###

#395,195,245

# Plot one ppx.
ppx=395
x <- 1:251
index <- x[ptsne.ppx==ppx]
# Find a value that makes the size function stable.
bdm.cost(big_list[[index]])
bdm.ptsne.plot(big_list[[index]], ptsne.cex = 2)
bdm.wtt.plot(big_list[[index]])
# The embedding function appears to reach stability the higher the ppx.
# However, for higher ppx, fewer clusters are formed, but the variability between threads increases.


# Get cluster labels
# get_ppx("allddi") 
# allddi,sen2in1,onlysen
get_ppx <- function(set_name){

  file <- file.path(paste0("data/1.processed/Broc2018_",set_name,"_set.csv"))
  dataset <- read_csv2(file) %>%
    filter(!is.na(ebw)&!is.na(ecr)&!is.na(seo)&!is.na(stm)&!is.na(pae)&!is.na(pau))
  
  big_list <- list_ppx[[set_name]]
  
  # Get cluster labels
  list_clusters <- list()
  for (i in 1:length(big_list)){
    list_clusters[[i]] <- bdm.labels(big_list[[i]])
  }
  number_of_clusters <- list_clusters |>
    map(unique) |> map_dbl(length)

  data.frame(index=1:251,
             ppx=ptsne.ppx,
             cluster_n=number_of_clusters)


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

  
write.csv(full_dataset,
          paste0("results/",set_name,"/","all_ppx_large_",set_name,".csv"),
          row.names = F)
}
