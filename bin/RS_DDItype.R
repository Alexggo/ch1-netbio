# Filter 2 groups to identify differences that are structural vs network based
# Set 1 SEN in all strains for both drugs (network differences). 
# Set 2 SEN/RES. Sensitive for both drugs in at least one species (structure and networks).
# Possible results
# Option 1. signal is the same, sensititvity doesn’t make any difference
# Option 2. Maybe effect is different.


library(tidyverse)

raw <- read_csv("data/1.processed/Broc2018_maindataset.csv")
x <-  raw |> 
  select(1:3,19:62)

y <- x |> 
  pivot_longer(names_to="RS_names",values_to = "RS_values",41:46) |> 
  pivot_longer(names_to="DDI_names",values_to = "DDI_values",16:21) |> 
  separate(RS_names,into=c("a","b","RS_names"),sep="_") |> 
  select(!c("a","b"))|> 
  separate(DDI_names,into=c("DDI_names","b"),sep="_") |> 
  select(!c("b")) |> 
  filter(RS_names==DDI_names)

y |> 
  mutate(cat=paste0(DDI_values,"-",RS_values)) |> 
  ggplot(aes(x=DDI_names,fill=as.factor(RS_values)))+
  theme_minimal()+
  geom_bar()+
  theme(legend.position = "bottom") +
  facet_wrap(~DDI_values)

y |> 
  mutate(cat=paste0(DDI_values,"-",RS_values)) |> 
  ggplot(aes(x=RS_names,fill=as.factor(DDI_values)))+
  theme_minimal()+
  geom_bar()+
  theme(legend.position = "bottom") +
  facet_wrap(~RS_values)

RR_non_add <- y |> 
  mutate(cat=paste0(DDI_values,"_",RS_values)) |> 
  filter(cat=="-1_+1"|cat=="-1_+1") |> pull(drug_pair)

RR_syn_res <- y |> 
  mutate(cat=paste0(DDI_values,"_",RS_values)) |> 
  filter(cat=="-1_1") |> pull(drug_pair)
# 62 DDI are synergistic and resistant

RR_ant_res <- y |> 
  mutate(cat=paste0(DDI_values,"_",RS_values)) |> 
  filter(cat=="1_1") |> pull(drug_pair)
# 11 DDI are antagonistic and resistant

RR_add_res <- y |> 
  mutate(cat=paste0(DDI_values,"_",RS_values)) |> 
  filter(cat=="0_1") |> pull(drug_pair)
#1927 are additive and resistant



SS_syn_sen <- y |> 
  mutate(cat=paste0(DDI_values,"_",RS_values)) |> 
  filter(cat=="-1_-1") |> pull(drug_pair)
# 289 DDI are synergistic and sensitive

SS_ant_sen <- y |> 
  mutate(cat=paste0(DDI_values,"_",RS_values)) |> 
  filter(cat=="1_-1") |> pull(drug_pair)
# 342 DDI are antagonistic and sensitive

SS_add_sen <- y |> 
  mutate(cat=paste0(DDI_values,"_",RS_values)) |> 
  filter(cat=="0_-1") |> pull(drug_pair)
#1165 are additive and sensitive





SR_syn_mix <- y |> 
  mutate(cat=paste0(DDI_values,"_",RS_values)) |> 
  filter(cat=="-1_0") |> pull(drug_pair)
# 300 DDI are synergistic and mixed

SR_ant_mix <- y |> 
  mutate(cat=paste0(DDI_values,"_",RS_values)) |> 
  filter(cat=="1_0") |> pull(drug_pair)
# 348 DDI are antagonistic and mixed

SR_add_mix <- y |> 
  mutate(cat=paste0(DDI_values,"_",RS_values)) |> 
  filter(cat=="0_0") |> pull(drug_pair)
#2541 are additive and mixed


conti <- matrix(c(289,300,62,
                  1165,2541,1927,
                  342,348,11),byrow = TRUE,nrow = 3)
colnames(conti) <- c("sen","hyb","res")
rownames(conti) <- c("syn","add","ant")



x_sq <- chisq.test(conti)

x_sq$observed < x_sq$expected
# ant and syn are less than expected when both targets are resistant
# add are less than expected in sen or sen/res targets.

#syn and ant are overrepresented in sen and hyb targets.
# add are overrepresented when both targets are resistant


sub <- x |> select(1:3,23:34)
# Set 1 SEN in all strains for both drugs (network differences). 
only_sen <- sub |> filter(SR_ebw_D1=="S"&SR_ecr_D1=="S"&SR_seo_D1=="S"&SR_stm_D1=="S"&SR_pae_D1=="S"&SR_pau_D1=="S") |> 
  filter(SR_ebw_D2=="S"&SR_ecr_D2=="S"&SR_seo_D2=="S"&SR_stm_D2=="S"&SR_pae_D2=="S"&SR_pau_D2=="S")
# 65 DDIs

# Set 2 SEN/RES. Sensitive for both drugs in at least one species (structure and networks).
set2 <- sub |> filter((SR_ebw_D1=="S"&SR_ebw_D2=="S")|(SR_ecr_D1=="S"&SR_ecr_D2=="S")|
                        (SR_seo_D1=="S"&SR_seo_D2=="S")|(SR_stm_D1=="S"&SR_stm_D2=="S")|
                        (SR_pae_D1=="S"&SR_pae_D2=="S")|(SR_pau_D1=="S"&SR_pau_D2=="S"))

set2.1 <- set2 |> filter(!drug_pair %in% only_sen$drug_pair)

only_sen_DDI <- only_sen |> select(drug_pair) |> pull()

set2.1_DDI <- set2.1 |> select(drug_pair) |> pull()

raw |> filter(drug_pair %in% only_sen_DDI) |> 
  write.csv("data/1.processed/Broc2018_sen_set.csv",row.names = FALSE)

raw |> filter(drug_pair %in% set2.1_DDI) |> 
  write.csv("data/1.processed/Broc2018_set2_set.csv",row.names = FALSE)
