#!/usr/bin/env Rscript

#' ---
#' title: "3.1.tSNE_BigMap"
#' author: "Alex Gil"
#' date: "1/26/2021"
#' output: html_document
#' 
#' 
#' 
## -------------------------------------------------------------------------
library(pacman)
p_load(bigMap,bigmemory,tidyverse,rlist)

args = commandArgs(trailingOnly=TRUE)

set.seed(1234)

file <- file.path("data/1.processed","Broc2018_set2_set.csv") #"Broc2018_maindataset.csv","Broc2018_sen_set.csv","Broc2018_set2_set.csv"
print(Sys.time())
#Remove DDI with at least one NA.
dataset <- read_csv(file) %>%
  filter(!is.na(ebw)&!is.na(ecr)&!is.na(seo)&!is.na(stm)&!is.na(pae)&!is.na(pau))

mat_red <- dataset %>% 
  select(ebw,ecr,stm,seo,pae,pau) %>%as.matrix

data <- mat_red

# input data big.matrix (a requisite to compute quality)
Dbm <- as.big.matrix(data, type = 'double')
# init model
m <- bdm.init('DDI', data)

# range of perplexities. 
# Large range for exploratory analysis.
ptsne.ppx <- seq(args[1], args[2], by = 10)
#Smaller range after plots are shown.
#ptsne.ppx <- seq(700, 900, by = 25)

# Should be around 1/3 of ptsne.ppx
pakde.ppx <- round(ptsne.ppx/3, 0)
# runs. 
# Number of threads should be as high as possible. Never higher than the number of hyper-threads.
# threads=cores*2
# Number of layers: the higher the most robust the result will be.

m.list1 <- lapply(seq_along(ptsne.ppx), function(p) {
  m.ppx <- bdm.ptsne(m, threads = 40, layers = 80, rounds = 9, whiten = 0, ppx = ptsne.ppx[p], info = 0)
  m.ppx <- bdm.pakde(m.ppx, threads = 40, ppx = pakde.ppx[p])
  m.ppx <- bdm.wtt(m.ppx)
  #m.ppx <- bdm.qlty(m.ppx, inp.data = Dbm, threads = 4, layers = 4, rounds = 2, ret.qlty = T, qm = 'kn', verbose = F)
  m.ppx
})

save(m.list1, file = paste0('results/set/tSNE_',args[1],"_",args[2],'.RData'))

print(Sys.time())
