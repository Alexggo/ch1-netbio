
#' ---
#' title: "1.1.Rearrangement"
#' author: "Alex Gil"
#
#Rearranges data from Brochado to generate maindataset.
## -----------------------------------------------------------------------------
library(pacman)
p_load(tidyverse,broom)

x <- read_csv(file.path("data/1.processed","2018_Broc_ED09C.csv")) %>% 
  separate(drug_pair,into=c("Drug1","Drug2"),sep="_") 
for (i in 1:dim(x)[1]){
  a <- x[i,]$Drug1
  b <- x[i,]$Drug2
  ab <- sort(c(a,b))
  x[i,]$Drug1 <- ab[1]
  x[i,]$Drug2 <- ab[2]
}
x <- x %>% mutate(drug_pair=paste0(Drug1,"_",Drug2))

y <- read_csv(file.path("data/1.processed","2018_Broc_SuppTable1.csv"))%>% select(1:5)
z <- read_csv(file.path("data/1.processed","2018_Broc_SuppTable2.csv")) %>%
  as_tibble() %>% 
  dplyr::rename(Drug1=drug1) %>% 
  dplyr::rename(Drug2=drug2) 

for (i in 1:dim(z)[1]){
  a <- z[i,]$Drug1
  b <- z[i,]$Drug2
  ab <- sort(c(a,b))
  z[i,]$Drug1 <- ab[1]
  z[i,]$Drug2 <- ab[2]
}

# Rearrange
y_red1 <- y %>% mutate(Drug1=drug,
                           drug_category1=drug_category,
                           targeted_cellular_process1=targeted_cellular_process,
                           use1=use,
                           code_3letter1=code_3letter) %>% 
  select(-c(1:5))
y_red2 <- y %>% mutate(Drug2=drug,
                           drug_category2=drug_category,
                           targeted_cellular_process2=targeted_cellular_process,
                           use2=use,
                           code_3letter2=code_3letter) %>% 
  select(-c(1:5))


y1 <- full_join(y_red2,x,by="Drug2")
y2 <- full_join(y_red1,y1,by="Drug1") 

z1 <- z %>% select(Drug1,Drug2,8:19) %>% 
  mutate(drug_pair=paste0(Drug1,"_",Drug2))


a <- full_join(y2,z1,by="drug_pair") %>% 
  select(-c(paste0("int_score_",c("ebw","ecr","seo","stm","pae","pau")))) %>% 
  filter(!is.na(drug_pair)) %>% 
  mutate(Drug1=Drug1.x,Drug2=Drug2.x) %>% 
  select(-c("Drug1.y","Drug2.y","Drug1.x","Drug2.x")) %>% 
  select(drug_pair,Drug1,Drug2,code_3letter1,code_3letter2,drug_category1,
         drug_category2,targeted_cellular_process1,targeted_cellular_process2,
         use1,use2,ebw,ecr,seo,stm,pae,pau,int_sign_ebw,int_sign_ecr,int_sign_seo,int_sign_stm,int_sign_pae,int_sign_pau) %>% 
  mutate(drug_category=ifelse(drug_category1==drug_category2,"Same","Different"),
         targeted_process=ifelse(targeted_cellular_process1==targeted_cellular_process2,"Same","Different"),
         use=ifelse(use1==use2,"Same","Different"))


d <- a %>% as.data.frame()
dd <- d %>% select(Drug1,Drug2)
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

# Add R/S info

res <- read_csv("data/0.rawdata/Broch_2018_278_MOESM3_ESM.csv") |> 
  select(1,6:23)

cols <- c("Drug","HDC_ebw_D","HDC_ecr_D",
          "HDC_stm_D","HDC_seo_D","HDC_pae_D","HDC_pau_D","LF_ebw_D",
          "LF_ecr_D","LF_stm_D",                 
          "LF_seo_D","LF_pae_D","LF_pau_D",                
          "SR_ebw_D","SR_ecr_D","SR_stm_D",                 
          "SR_seo_D","SR_pae_D","SR_pau_D")

colnames(res) <- cols

res1 <- res
colnames(res1) <- paste0(cols,"1")

res2 <- res
colnames(res2) <- paste0(cols,"2")

e1 <- left_join(e,res1,by="Drug1")
e2 <- left_join(e1,res2,by="Drug2")

e3 <- e2 |> 
  mutate(SR_ebw=paste0(SR_ebw_D1,SR_ebw_D2),
         SR_ecr=paste0(SR_ecr_D1,SR_ecr_D2),
         SR_stm=paste0(SR_stm_D1,SR_stm_D2),
         SR_seo=paste0(SR_seo_D1,SR_seo_D2),
         SR_pae=paste0(SR_pae_D1,SR_pae_D2),
         SR_pau=paste0(SR_pau_D1,SR_pau_D2)) |> 
  mutate(SR_score_ebw=ifelse(SR_ebw=="RR",1,
                             ifelse(SR_ebw=="SR"|SR_ebw=="RS",0,
                                    ifelse(SR_ebw=="SS",-1,"NA"))))|> 
  mutate(SR_score_ecr=ifelse(SR_ecr=="RR",1,
                             ifelse(SR_ecr=="SR"|SR_ecr=="RS",0,
                                    ifelse(SR_ecr=="SS",-1,"NA"))))|> 
  mutate(SR_score_stm=ifelse(SR_stm=="RR",1,
                             ifelse(SR_stm=="SR"|SR_stm=="RS",0,
                                    ifelse(SR_stm=="SS",-1,"NA"))))|> 
  mutate(SR_score_seo=ifelse(SR_seo=="RR",1,
                             ifelse(SR_seo=="SR"|SR_seo=="RS",0,
                                    ifelse(SR_seo=="SS",-1,"NA"))))|> 
  mutate(SR_score_pae=ifelse(SR_pae=="RR",1,
                             ifelse(SR_pae=="SR"|SR_pae=="RS",0,
                                    ifelse(SR_pae=="SS",-1,"NA"))))|> 
  mutate(SR_score_pau=ifelse(SR_pau=="RR",1,
                             ifelse(SR_pau=="SR"|SR_pau=="RS",0,
                                    ifelse(SR_pau=="SS",-1,"NA")))) |> 
  mutate(SR_score_ebw=as.numeric(SR_score_ebw)) |> 
  mutate(SR_score_ecr=as.numeric(SR_score_ecr)) |> 
  mutate(SR_score_stm=as.numeric(SR_score_stm)) |> 
  mutate(SR_score_seo=as.numeric(SR_score_seo)) |> 
  mutate(SR_score_pae=as.numeric(SR_score_pae)) |> 
  mutate(SR_score_pau=as.numeric(SR_score_pau)) |> 
  rowwise() |> 
  mutate(SR_score_total=SR_score_ebw+SR_score_ecr+SR_score_seo+SR_score_stm+
                            SR_score_pae+SR_score_pau)|>
  mutate(drug_cat = ifelse(drug_category == "Same", 1, 0)) |>
  mutate(targ = ifelse(targeted_process == "Same", 1, 0)) |>
  mutate(us = ifelse(use == "Same", 1, 0)) |>
  mutate(ebw_g = ifelse(int_sign_ebw == "Antagonism",  +1,
                        ifelse(int_sign_ebw == "Additivity", 0,
                               ifelse(int_sign_ebw == "Synergy", -1, NA)))) |>
  mutate(ecr_g = ifelse(int_sign_ecr == "Antagonism",  +1,
                        ifelse(int_sign_ecr == "Additivity", 0,
                               ifelse(int_sign_ecr == "Synergy", -1, NA)))) |>
  mutate(seo_g = ifelse(int_sign_seo == "Antagonism",  +1,
                        ifelse(int_sign_seo == "Additivity", 0,
                               ifelse(int_sign_seo == "Synergy", -1, NA)))) |>
  mutate(stm_g = ifelse(int_sign_stm == "Antagonism",  +1,
                        ifelse(int_sign_stm == "Additivity", 0,
                               ifelse(int_sign_stm == "Synergy", -1, NA)))) |>
  mutate(pae_g = ifelse(int_sign_pae == "Antagonism",  +1,
                        ifelse(int_sign_pae == "Additivity", 0,
                               ifelse(int_sign_pae == "Synergy", -1, NA)))) |>
  mutate(pau_g = ifelse(int_sign_pau == "Antagonism",  +1,
                        ifelse(int_sign_pau == "Additivity", 0,
                               ifelse(int_sign_pau == "Synergy", -1, NA)))) |>
  rowwise() |>
  mutate(sum_g = ebw_g + ecr_g + seo_g + stm_g + pae_g + pau_g) |> 
  mutate(only_res=ifelse(SR_score_total==6,TRUE,FALSE)) |> 
  mutate(only_sen=ifelse(SR_score_total==-6,TRUE,FALSE)) |> 
  mutate(only_res_or_sen=ifelse(only_res==TRUE|only_sen==TRUE,TRUE,FALSE))

e4 <- e3 |> select(drug_pair,Drug1,Drug2,drugdrug,
                   code_3letter1,code_3letter2,
                   drug_cat,drug_category1,drug_category2,categorycategory,
                   targ,targeted_cellular_process1,targeted_cellular_process2,processprocess,
                   us,use1,use2,useuse,
                   12:23,83:89,
                   43:48,61:79,
                   31:60)

write_csv2(e4,"data/1.processed/Broc2018_allddi_set.csv")
