#' ---
#' title: "1.CIDlookup.Rmd"
#' author: "Alex Gil"
#' date: "1/25/2021"
#' output: html_document
#' ---
#' 
#' # 1.CIDlookup.Rmd
#' 
#' Script that obtains the CID values for a given list of chemical compounds, such as Spitzer 2015.
#' 
## -----------------------------------------------------------------------------
library("plyr")
library("tidyverse")
library("webchem")



#'
#' 
## -----------------------------------------------------------------------------
filepath <- "data/0.rawdata"
filename <- "2018_Broc_SuppTable2.csv" 

input <- file.path(filepath,filename)

tab_compound <- read_csv(input)
tab_compound <- tab_compound[-1,]
colnames(tab_compound)[1:2] <- c("drug1","drug2")
names <- c(tab_compound$drug1,tab_compound$drug2) %>% 
  unique() %>% 
  sort()  %>% as.data.frame()
colnames(names) <- "names"

names <- names %>% as_tibble %>% mutate(cid = map(names,get_cid))


final_table <- names  %>% 
  mutate(cid = map(cid,select,2)) %>%
  mutate(cid=map(cid,head,1)) %>% 
  mutate(cid=map(cid,as.numeric)) %>% 
  unnest(cid) %>% 
  as.data.frame()

write.csv(final_table,"data/3.Targets_NetworkDistance/DrugTargets0_pubchemid.csv",row.names = FALSE)



