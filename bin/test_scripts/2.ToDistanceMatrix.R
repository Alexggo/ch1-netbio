#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#' ---
#' title: "2.ToDistanceMatrix"
#' author: "Alex Gil"
#' date: "1/25/2021"
#' output: html_document
#' ---
#' # 2.ToDistanceMatrix.Rmd
#' 
#' Takes a table a converts it to a distance matrix
#' 
## -----------------------------------------------------------------------------


library("plyr")
library("tidyverse")

file <- args[1]
tab <- read_csv(file=file.path("data/processed","2018_Broc_ED09C.csv")) #file


tabnew <- tab %>% 
  select(-1)  %>%
  t() %>% 
  dist()%>% 
  as.matrix() %>% 
  as.data.frame()

tabnew <- tabnew[ , order(names(tabnew))]
tabnew$Sp1 <- rownames(tabnew)

tabnew <- tabnew %>%
  arrange(Sp1) %>% 
  select(-Sp1)

tabnew[lower.tri(tabnew)] <- NA

tabnew$Sp1 <- colnames(tabnew)

tabnew1 <- tabnew %>%
  select(dim(tabnew)[2],1:(dim(tabnew)[2]-1))

tabnew2 <- tabnew1  %>% 
  pivot_longer(names_to="Sp2",values_to="Chemdist",2:(dim(tabnew)[2]))%>%
  filter(!is.na(Chemdist)) %>% 
  filter(Chemdist!=0)

write.csv(tabnew2,file.path("data/processed",paste0(file,".mat.csv")),row.names = F)

