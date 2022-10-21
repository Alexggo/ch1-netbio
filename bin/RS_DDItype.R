library(tidyverse)

x <- read_csv("data/1.processed/Broc2018_maindataset.csv") |> 
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
