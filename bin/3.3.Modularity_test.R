#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#Rscript --vanilla 3.3.Modularity_test.R set_name ppx ppx


#' title: "3.2.Modularity_test"
#' author: "Alex Gil"
#' date: "1/26/2021"
#' output: html_document

# Run with:
# Rscript bin/3.3.Modularity_test.R allddi ppx_205 ppx_1455

library(pacman)
p_load(geomorph,broom,ape,ggx,furrr,tictoc,dplyr,tidyr,readr,purrr)


availableCores()
plan(multisession, workers = 28)

# Read phylogenetic tree
species <- c("Escherichia_coli_K-12_ebw",
             "Escherichia_coli_O8_IAI1_ecr",
             "Salmonella_enterica_serovar_Typhimurium_LT2_stm",
             "Salmonella_enterica_serovar_Typhimurium_14028S_seo",
             "Pseudomonas_aeruginosa_PAO1_pae",
             "Pseudomonas_aeruginosa_UCBPP-PA14_pau")

treefile <- file.path('data/2.2.Phylogenetics_Bayesian',"Phylo-p1.nwk.tree")
tree_nw <- read.nexus(treefile)


# ALL DDI. Good tsnes between 205-1455.
# SEN2IN1. Good tsnes between 115-455.
# ONLYSEN. All tsnes are bad due to the small number of DDIs.


#get_ppx("allddi","ppx_205","ppx_1455")
#get_ppx("sen2in1","ppx_115","ppx_455")

get_ppx <- function(set_name,start,end){
  print(paste(start,'to',end))
full_dataset <- read_csv(paste0("results/",set_name,"/","all_ppx_large_",set_name,".csv"))

mat <- full_dataset[,19:24]%>% as.matrix() |> t()

## Testing modularity of clustering solutions.
mod_test <- list()
labels <- list()
cols <- 1:dim(full_dataset)[2]
sta_c <- cols[colnames(full_dataset)==start]
end_c <- cols[colnames(full_dataset)==end]

for (i in sta_c:end_c){
  labels[[i-sta_c+1]] <- full_dataset[,i] |> pull()
}

number_of_clusters <- labels |> map(unique) |> map_dbl(length)

test_mod <- function(g){
  mod <- phylo.modularity(mat,partition.gp = g,
                   phy=tree_nw,print.progress=T,
                   CI = TRUE)
}

tic()
mod_test <- labels |>
  lapply(test_mod)
toc()

file_name <- paste0('results/',set_name,'/modularity_tests_',
                    set_name,start,"-",end,'.RData')
print(file_name)
mod_test |> save(file = file_name)
}

get_ppx(args[1],args[2],args[3])

