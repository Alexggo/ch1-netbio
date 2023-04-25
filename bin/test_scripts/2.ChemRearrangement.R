#!/usr/bin/env Rscript


#' ---
#' title: "2.ChemRearrangements"
#' author: "Alex Gil"
#' date: "1/25/2021"
#' output: html_document
#' ---
#' # 1.ChemRearrangement.Rmd
#' 
#' Script that rearranges the chemical results from different papers into distance matrices.
#' 
## -----------------------------------------------------------------------------
library("plyr")
library("tidyverse")

#' 
#' ## Spitzer 2011: Fungi
#' 
#' ### Contents:
#' 
#'   * DataHetES
#'   * DataSetTransform
#'   * Z-scores 
#'   * S1. High-throughput screen data. Spit2011S1
#'   * S2. Minimum inhibitory concentrations of hit compounds
#'   * S3. Fractional inhibitory concentration index of hit compounds
#'   * S4. P-values for mean Z-scores of fluconazole-specific deletion strains in indicated chemical-genetic profiles.
#'   * S5. Z-scores for haploid deletion strains that are significant in at least one drug screen
#'   * S6. P-values for the CGS and PPP simulations for six drug profiles
#'   * S7. MIC and FICI data for the fluconazole - sertraline drug combination against drug resistant Candida strains
#'   * FigS1. Normalization of residual growth data
#'   * FigS2. Activity of hit compounds in Prestwick library screens
#'   * FigS3. Chemical structures and application classification of 17 compounds selected for synergy confirmation experiments
#'   * FigS4. Assessment of primary high throughput screen performance
#'   * FigS5. Confirmation of hits from haplo-insufficiency chemical-genetic profiles. Average Z-score
#'   * FigS6. Global heatmap representing all haploid deletion strain chemical-genetic profiles in this study
#'   * FigS7. Sensitivity of all haploid deletion strains treated with each of the syncretic synergizers
#'   * FigS8.  Effects of syncretic compounds on membrane integrity
#'   * FigS9. Significance of genetic interaction enrichments between fluconazole and synergistic drugs by parallel pathway permutation test
#' 
#' Input the files.
#' 
## -----------------------------------------------------------------------------
#Spit2011_S1
df2011 <- read_csv("data/rawdata/2011_Spit_S1.csv")
colnames(df2011) <- c("species","plate","column","row","compound","prestwick","flc_pre","inh",
                      "residuals","hit")

df2011$species <- revalue(df2011$species, c("S.cerevisiae"="sce",
                                            "C.albicans"="cal",
                                            "C.neoformans"="cne",
                                            "C.gattii"="cgi"))

df2011.1 <- df2011 %>% 
  select(species,compound,inh) %>% 
  pivot_wider(names_from = species,values_from=inh) %>% 
  mutate(cal=as.numeric(cal),
         cgi=as.numeric(cgi),
         cne=as.numeric(cne),
         sce=as.numeric(sce))

write.csv(file="data/processed/2011_Spit_S1.1.csv",df2011.1,row.names=FALSE)

#Spit2011_S3. #FICI values with respect to fluconazole
FICI <- read.csv(file = "data/rawdata/2011_Spit_S3_FICI.csv",header = T)
colnames(FICI) <- c("compound","sce","cal","cne","cgi")

write.csv(file="data/processed/2011_Spit_S3.1.csv",FICI,row.names = F)


#' ## Spitzer 2015: Fungi
#' 
#' ### Contents:
#' 
#'  * S2 GrowthMat Script 2015
#'  		-   Rscript Spit2015E_S2GrowthMat.R
#'  		-   Rscript Spit2015E_S2GrowthMatCID.R
#'  * S3. Table Antifungal+Compound
#'  * S5. FICImat script
#'  		-   less Spit2015S5.csv | Rscript Spit2015E.S5CalcFICImat.R Spit2015ES5distFICI.txt
#'  * S6. Homozygous deletion profiles and haploinsufficiency profile. Z scores + 3 replicates
#'  * S8. Primer names and sequences
#' 
#' Input the files.
## -----------------------------------------------------------------------------
##2015_Spit_S2
df2015 <- read_csv("data/rawdata/2015_Spit_S2.csv") %>% 
  select(-c(2:4))
colnames(df2015) <- c("species",'compound','DMSO','AmpB','Ben',
                      'CF','Cyp','FL','TB')

df2015$species <- revalue(df2015$species,c("S.cerevisiae"="sce",
                                           "C.albicans"="cal",
                                           "C.neoformans"="cne",
                                           "C.gattii"="cgi",
                                           "S.pombe"="spo"))

df2015 <- df2015 %>% pivot_longer(names_to="Treatment",values_to = "Growth",3:9)%>% 
  arrange(compound) %>% 
  filter(!is.na(compound)) %>% 
  mutate(compound=tolower(as.character(compound)))

df2015.1 <- df2015 %>% 
  mutate(DrugPair = paste0(compound,"_",Treatment)) %>% 
  select(1,5,4) %>% 
  rename(compound=DrugPair) %>% 
  group_by(species,compound)%>%
  summarise(Growth=mean(Growth)) %>% 
  pivot_wider(names_from=compound,values_from=Growth)

df2015.1 <- df2015.1[ , apply(df2015.1, 2, function(x) !any(is.na(x)))]

df2015.1 <- df2015.1 %>% as.data.frame()
rownames(df2015.1) <- df2015.1[,1]

df2015.1 <- df2015.1 %>% select(-1)

df2015.1 <- df2015.1 %>% t() %>% as.data.frame()

df2015.1$compound <- rownames(df2015.1)

df2015.1 <- df2015.1 %>% select(5,1:4)

rownames(df2015.1) <- NULL

write.csv(df2015.1,"data/processed/2015_Spit_S2.1.csv",row.names = FALSE)


# Spit2015.S5. ACM compound chemical interactions with antifungals. 
# Calculated Fractional Inhibitory Concentration Index. Related to Figure 3.
df2015.5 <- read.csv("data/rawdata/2015_Spit_S5.csv", sep=",")
colnames(df2015.5) <- c("species","compound","antifungal","cid","MIC_ACM","FIC_ACM","MIC_Antifungal","FIC_Antifungal","FIC_index")

df2015.5$species <- revalue(df2015.5$species, c("S.cerevisiae"="sce",
                                              "C.albicans"="cal",
                                              "C.neoformans"="cne",
                                              "C.gattii"="cgi",
                                              "S.pombe"="spo"))
#Select FIC columns
df2015.5 <- df2015.5 %>% mutate(drugcomb=paste0(antifungal,"_",compound)) %>% 
  select(1,6:10)%>%
  group_by(species,drugcomb)%>%
  summarise(FIC_ACM =mean(FIC_ACM),
            FIC_Antifungal=mean(FIC_Antifungal),
            FIC_index=mean(FIC_index)) %>% 
  select(1,2,5)%>% 
  pivot_wider(names_from=drugcomb,values_from=FIC_index)%>% 
  as.data.frame()

rownames(df2015.5) <- df2015.5[,1]
df2015.5 <-  df2015.5%>%
  select(-1) %>% t() %>% 
  as.data.frame()

df2015.5$compound <- rownames(df2015.5)
df2015.5 <- df2015.5 %>% select(5,1:4)

rownames(df2015.5) <- NULL

write.csv(file="data/processed/2015_Spit_S5.1.csv",df2015.5,row.names=FALSE)


#' 
#' 
#' ## Brochado 2018: Bacteria
#' 
#'   * ED2
#'     - Panel B (Growth of single query drug, growth double drug query+spectominomycin)
#'     - Panel C (Bliss scores, medians and interquantile range)
#'  			+ Rscript Bro2018ED2Script.R
#'   * ED3
#'     - Panel A (Pearson Correlations of single drugs and double drug)
#'     - Panel B (Pearson Coefficient of correlation)
#'     - Panel C (Nothing: refer to supplementary file 1)
#'     - Panel D (DrugPairs: Expected fitness and Bliss scores)
#'     - Panel E (Lowest single fitness, nr interactions per drug and strain)
#'     - Panel F (nr interactions per drug)
#'   * ED4
#'     - Panel A (screen and benchmarking)
#'     - Panel B (Benchmarking combinations)
#'     - Panel C (interactions)
#'     - Panel D (RP rate, FP rate)
#'     - Panel E (Fitness for different doses of drugpairs in Salmonella and E.Coli)
#'   * ED5
#'     - Panel A (Benchmarking and nr interactions)
#'  			+ PAO1 (Bliss score and expected fitness)
#'  			+ PA14 (Bliss score and expected fitness)
#'  			+ PAO-replicate 1+2 (checkerboard values)
#'  			+ PA12-replicate 1+2 (checkerboard values)
#'   * ED6
#'     - Panel A (Benchmarking and nr interactions)
#'  			+ SL-LT2 (Bliss score and expected fitness)
#'  			+ SL-14028S (Bliss score and expected fitness)
#'  			+ SL-LT2-replicate 1+2 (checkerboard values)
#'  			+ SL-14028S-replicate 1+2 (checkerboard values)
#'   * ED7
#'     - Panel A (categories of drugs 1 and 2 and number of interactions)
#'     - Panel B (categories of drugs 1 and 2 and number of interactions)
#'     - Panel C (categories of drugs 1 and 2 and number of interactions)
#'     - Panel D (categories of drugs 1 and 2 and number of interactions)
#'     - Panel E (Drug general target)
#'   * ED8
#'     - Panel B (Intracellular concentration ratio)
#'     - Panel C (Interaction Scores)
#'     - Panel D (Counts and Drugs)
#'   * ED9
#'     - Panel A (interaction scores strains comparisons ST)
#'     - Panel B (interaction scores strains comparisons PA)
#'     - Panel C (interaction scores for drug pairs each bliss score)
#'     - Panel D (Interaction and Conservation)
#'     - Panel E (Monochromaticity index)
#'     - Panel F (Monochromaticity index per class)
#'   * ED10
#'     - Sheet (median of bliss interaction distribution, drug pairs per strain)
#'  			+ Rscript Bro2018ED10Script.R
#'   * ED11
#'     - Panel A (Please refer to supplementary file 8 for checkerboards)
#'     - Panel B (Survival and standard deviation, time hpi)
#'   * ED12
#'     - Panel A (Streptomycin MIC)
#'     - Panel B (MICs per strain)
#'     - Panel C (refer to supp file 7)
#'     - Panel D (Streptomyvin MIC)
#'     - Panel E (Growth after 8h with different dosis of streptomycin)
#'     - Panel F (Fold change and treatment)
#'   * Fig1
#'     - Panel A (synergy and antagonism, detected, interactions)
#'     - Panel B (category drug1-drug2, nr interactions, direction)
#'     - Panel C&E (Drug category, drug general target)
#'     - Panel D (General target drug1:drug2, nr interactions, direction)
#'   * Fig2
#'     - Panel A (drug pair, interaction scores E.coli)
#'     - Panel B (nr interactions, conserved...)
#'     - Panel C (Drug interactions, species conservation)
#'     - Panel D (Interaction, conserved/nonconserved)
#'   * Fig3
#'     - Panel A (refer to ED10)
#'     - Panel B (Strain, AcrA levels, treatment)
#'     - Panel C (Strain, marA expre, treatment)
#'     - Panel D (Strains, MICs, treatment)
#'     - Panel E (Strains, MICs, treatment)
#'   * Fig4
#'     - Panel A (refer to file8 for in vitro checkerboards)
#'     - Panel B (%survival, standard deviation, tiem hpi, treatment)
#' 
#' Input the files.
## -----------------------------------------------------------------------------
# 2018_Broc_ED02.csv
tab <- read.csv("data/rawdata/2018_Broc_ED02.csv",header = F)

# 2018_Broc_ED02.csv
for(i in 1:19){
  colnames(tab)[i] <- paste(tab[1,i],"-",tab[2,i])
}

tab <- tab[-c(1,2),]
colnames(tab)[1] <- "Drug Pair"
tab <- tab %>% pivot_longer(2:19,names_to="Variable",values_to="Value")

tab <- tab %>% separate(col=Variable,into=c("spst","stat"),sep="-") %>% 
  mutate(Value=as.numeric(Value))

tab$spst <- revalue(tab$spst, c("PAO1 "="pae",
                                "PA14 "="pau",
                                "ST LT2 "="stm",
                                "ST 14028s "="seo",
                                "EC BW "="ebw",
                                "EC iAi1 "="ecr"))

tab$stat <- revalue(tab$stat, c(" Median Q2"="Q2",
                                " Q1"="Q1",
                                " Q3"="Q3"))

write.csv(tab,"data/processed/2018_Broc_ED02.1.stats.csv",row.names = F)

# Just the medians.
m1 <- tab %>% 
  filter(stat=="Q2") %>% 
  select(1,2,4) %>% 
  pivot_wider(names_from = spst,values_from=Value) %>% 
  rename(compound="Drug Pair") %>% 
  filter(!is.na(pae)) %>% 
  filter(!is.na(pau))

write.csv(m1,"data/processed/2018_Broc_ED02.1.csv",row.names = F)

# Calculate Stats.
## Remove NAs. Mean/Variance estimations.
#Estimated mean of the sample from Wan et al. (2014).
# SD estimation. Cochrane Handbook (2008).
table.stats <- read.csv("data/processed/2018_Broc_ED02.1.stats.csv")

table1 <- table.stats %>% 
  filter(!is.na(Value)) %>% 
  pivot_wider(names_from = stat,values_from=Value) %>% 
  mutate(mean=(Q1 + Q2 + Q3)/3,
         sd = (Q3 - Q1)/1.35)

write.csv(table1,"data/processed/2018_Broc_ED02.2.stats.csv",row.names = F)


#number of interactions matrix. Broc_ED03 PanelF.
Bro2018.S3.1 <- read.csv(file="data/rawdata/2018_Broc_ED03F.nr.csv")
colnames(Bro2018.S3.1)[1] <- "compound"

Bro2018.S3.1 <- Bro2018.S3.1 %>%
  pivot_longer(2:7,names_to="spst",values_to="nr")

Bro2018.S3.1$spst <- revalue(Bro2018.S3.1$spst, c("PAO1"="pae",
                                "PA14"="pau",
                                "ST.LT2"="stm",
                                "ST.14028s"="seo",
                                "Ecoli.BW"="ebw",
                                "Ecoli.iAi1"="ecr"))

Bro2018.S3.1 <- Bro2018.S3.1 %>% 
  pivot_wider(names_from="spst",values_from="nr")

write.csv(Bro2018.S3.1,"data/processed/2018_Broc_ED03.2.csv",row.names = F)

