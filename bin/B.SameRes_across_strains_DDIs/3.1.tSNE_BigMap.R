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

file <- file.path("data/1.processed",
                  "Broc2018_maindataset.csv")
print(Sys.time())
#Remove DDI with at least one NA.
dataset <- read_csv(file) %>%
  filter(!is.na(ebw)&!is.na(ecr)&!is.na(seo)&!is.na(stm)&!is.na(pae)&!is.na(pau))

#Filter only drugs with the same type of resistance level across strains. Sup Tab1.
# All R, and all S.
same_res <- c("Imipenem","Aztreonam","Cephalexin","Teicoplanin",
             "Bacitracin","Metronidazole","Clofazimine",
             "Doxorubicin","Mitomycin C","Ciprofloxacin","Levofloxacin",
             "Sulfamonomethoxine","Cerulenin","Nisin","Daptomycin",
             "Nonactin","Paraquat","Linezolid","Tigecycline",
             "Spiramycin","Erythromycin","Clarithromycin",
             "Caffeine","Reserpine","Procaine","Theophylline",
             "Verapamil","Meropenem","Triclosan",
             "Gentamicin","Berberine","Metformin","Diclofenac")

dataset <- dataset |> 
  filter(Drug1 %in% same_res) |> 
  filter(Drug2 %in% same_res)

mat_red <- dataset %>% 
  select(ebw,ecr,stm,seo,pae,pau) %>%as.matrix

data <- mat_red

# input data big.matrix (a requisite to compute quality)
Dbm <- as.big.matrix(data, type = 'double')
# init model
m <- bdm.init('DDI', data)

# range of perplexities. 
# Large range for exploratory analysis.
ptsne.ppx <- c(seq(5,50,10),seq(50, 2120, by = 82))
#Smaller range after plots are shown.
#ptsne.ppx <- seq(700, 900, by = 25)

# Should be around 1/3 of ptsne.ppx
pakde.ppx <- round(ptsne.ppx/3, 0)
# runs. 
# Number of threads should be as high as possible. Never higher than the number of hyper-threads.
# threads=cores*2
# Number of layers: the higher the most robust the result will be.

m.list1 <- lapply(seq_along(ptsne.ppx), function(p) {
  m.ppx <- bdm.ptsne(m, threads = 80, layers = 80, rounds = 9, whiten = 0, ppx = ptsne.ppx[p], info = 0)
  m.ppx <- bdm.pakde(m.ppx, threads = 40, ppx = pakde.ppx[p])
  m.ppx <- bdm.wtt(m.ppx)
  #m.ppx <- bdm.qlty(m.ppx, inp.data = Dbm, threads = 4, layers = 4, rounds = 2, ret.qlty = T, qm = 'kn', verbose = F)
  m.ppx
})

save(m.list1, file = 'results/SameRes_across_strains_DDIs/RObjects/mlist1.RData')

print(Sys.time())

