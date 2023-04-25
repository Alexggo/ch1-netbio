#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)



#' ---
#' title: "1.PhyloRearrangement.Rmd"
#' author: "Alex Gil"
#' date: "1/25/2021"
#' output: html_document
#' ---
#' 
#' # 1.PhyloRearrangement.Rmd
#' 
#' Script that rearranges the phylogenetic results from PhySpeTree into distance matrices.
#' 
## -----------------------------------------------------------------------------
##HOW TO RUN: Rscript 1.PhyloRearrangement.R pathtofile

library("plyr")
library("tidyverse")

#' 
#' # Input your data in the file path, and run the chunk.
#' 
#' ../data/PhySpeTree/bac_species.txt.hcp
#' ../data/PhySpeTree/bac_species.txt.srna
#' ../data/PhySpeTree/fun_species.txt.hcp
#' ../data/PhySpeTree/fun_species.txt.srna
#' 
## -----------------------------------------------------------------------------




filepath <- args[1]
taxa <- str_extract(filepath,".../*_")
seq_type <- str_extract(filepath,regex("...$"))

input <- file.path(filepath,"iqtree.tree.mldist")

matdist <- read_delim(input,skip=1,delim=" ",col_names=F)
matdist <- matdist[, colSums(is.na(matdist)) != nrow(matdist)]
matdist <- matdist %>% as.data.frame()
colnames(matdist) <- c("Sp1",as.character(matdist[,1]))
matdist <- matdist[ , order(names(matdist))]
rownames(matdist) <- matdist$Sp1

matdist <- matdist %>%
  arrange(Sp1)
matdist <- matdist %>%
  select(-Sp1)

matdist[lower.tri(matdist)] <- NA
matdist$Sp1 <- rownames(matdist)
d <- dim(matdist)
matdist <- matdist %>%
  select(d[2],1:(d[2]-1))

matdist1 <- matdist %>%
  as_tibble()  %>%
  gather(key="Sp2",value="PhyloDist",2:(dim(matdist)[2]))%>%
  filter(!is.na(PhyloDist)) %>% 
  mutate(PhyloDist=as.numeric(PhyloDist))

output <- file.path("../data/PhySpeTree",paste0("dist_",taxa,seq_type,".csv"))
write.csv(file=output,matdist1,row.names=FALSE)

