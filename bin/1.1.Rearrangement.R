
#' ---
#' title: "1.1.Rearrangement"
#' author: "Alex Gil"
#
#Rearranges data from Brochado to generate maindataset.
## -----------------------------------------------------------------------------
library("plyr")
library("tidyverse")
library("broom")

x <- read_csv(file.path("data/1.processed","2018_Broc_ED09C.csv"))
y <- read_csv(file.path("data/1.processed","2018_Broc_SuppTable1.csv"))%>% select(1:5)
z <- read_csv(file.path("data/1.processed","2018_Broc_SuppTable2.csv"))

# Rearrange
y_red1 <- y %>% mutate(drug1=drug,
                           drug_category1=drug_category,
                           targeted_cellular_process1=targeted_cellular_process,
                           use1=use,
                           code_3letter1=code_3letter) %>% 
  select(-c(1:5))
y_red2 <- y %>% mutate(drug2=drug,
                           drug_category2=drug_category,
                           targeted_cellular_process2=targeted_cellular_process,
                           use2=use,
                           code_3letter2=code_3letter) %>% 
  select(-c(1:5))


x <- x %>% 
  separate(drug_pair, c("drug1", "drug2"), "_") %>% 
  mutate(drug_pair=paste0(drug1,"_",drug2))

y1 <- full_join(y_red2,x,by="drug2")
y2 <- full_join(y_red1,y1,by="drug1") 

z1 <- z %>% select(drug1,drug2,8:19) %>% 
  mutate(drug_pair=paste0(drug1,"_",drug2))


a <- full_join(y2,z1,by="drug_pair") %>% 
  select(-c(paste0("int_score_",c("ebw","ecr","seo","stm","pae","pau")))) %>% 
  filter(!is.na(drug_pair)) %>% 
  mutate(drug1=drug1.x,drug2=drug2.x) %>% 
  select(-c("drug1.y","drug2.y","drug1.x","drug2.x")) %>% 
  select(drug_pair,drug1,drug2,code_3letter1,code_3letter2,drug_category1,
         drug_category2,targeted_cellular_process1,targeted_cellular_process2,
         use1,use2,ebw,ecr,seo,stm,pae,pau,int_sign_ebw,int_sign_ecr,int_sign_seo,int_sign_stm,int_sign_pae,int_sign_pau) %>% 
  mutate(drug_category=ifelse(drug_category1==drug_category2,"Same","Different"),
         targeted_process=ifelse(targeted_cellular_process1==targeted_cellular_process2,"Same","Different"),
         use=ifelse(use1==use2,"Same","Different"))




d <- a %>% as.data.frame()
dd <- d %>% select(drug1,drug2)
cc <- d %>% select(drug_category1,drug_category2)
pp <- d %>% select(targeted_cellular_process1,targeted_cellular_process2)
uu <- d %>% select(use1,use2)

for (i in 1:dim(d)[1]){
  drugdrug <- paste0(sort(dd[i,])[1],"-",sort(dd[i,])[2])
  d$drugdrug[i] <- drugdrug
  categorycategory <- paste0(sort(cc[i,])[1],"-",sort(cc[i,])[2])
  d$categorycategory[i] <- categorycategory
  processprocess <- paste0(sort(pp[i,])[1],"-",sort(pp[i,])[2])
  d$processprocess[i] <- processprocess
  useuse <- paste0(sort(uu[i,])[1],"-",sort(uu[i,])[2])
  d$useuse[i] <- useuse

}

e <- d %>% mutate(int_sign_ebw=ifelse(is.na(int_sign_ebw),"Additivity",int_sign_ebw)) %>% 
  mutate(int_sign_ecr=ifelse(is.na(int_sign_ecr),"Additivity",int_sign_ecr)) %>% 
  mutate(int_sign_seo=ifelse(is.na(int_sign_seo),"Additivity",int_sign_seo)) %>% 
  mutate(int_sign_stm=ifelse(is.na(int_sign_stm),"Additivity",int_sign_stm)) %>% 
  mutate(int_sign_pae=ifelse(is.na(int_sign_pae),"Additivity",int_sign_pae)) %>% 
  mutate(int_sign_pau=ifelse(is.na(int_sign_pau),"Additivity",int_sign_pau)) 

write.csv(e,"data/1.processed/Broc2018_maindataset.csv",row.names = F)

