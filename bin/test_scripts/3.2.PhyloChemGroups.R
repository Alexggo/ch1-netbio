#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


#' ---
#' title: "3.2.PhyloChemGroups.Rmd"
#' author: "Alex Gil"
#' date: "1/26/2021"
#' output: html_document
#' ---
#' 
#' # 3.PhyloChem.Rmd
#' 
#' Script that compares phylogenetic and chemical distances between species of fungi and bacteria.
#' HOW TO USE? ARGUMENT 1 is the phylogenetic distance matrix, ARGUMENT 2 is the chemical distance matrix.
## ----include=FALSE------------------------------------------------------------
library("plyr")
library("tidyverse")
library("broom")
library("ggrepel")
library("hrbrthemes")
library("ggforce")
library("scales")
library("plotly")

#' 
#' ## Input phylogenetic distance matrices. args[1]
#'   + Bacteria
#'     - dist_bac_hcp.csv
#'     - dist_bac_rna.csv
#'   + Fungi
#'     - dist_fun_hcp.csv
#'     - dist_fun_rna.csv
#' 
#' ## Input chemical distance matrix. args[2]
#'   + Bacteria
#'     - Bliss
#'       2018_Broc_ED02.1.csv.mat.csv
#'     - NR
#'       2018_Broc_ED03.2.csv.mat.csv
#' 
#'   + Fungi
#'     - Sp2011.s1 (flc_prw growth)
#'       2011_Spit_S1.1.csv.mat.csv
#'     - Sp2011.s3 FICI(17 x fluconazole)
#'       2011_Spit_S3.1.csv.mat.csv
#'     - Sp2015.s2 Growth(17 x 8 antifungals)
#'       2015_Spit_S2.1.csv.mat.csv
#'     - Sp2015.s5 FICI(17 x 2 antifungals)
#'       2015_Spit_S5.1.csv.mat.csv  
## -----------------------------------------------------------------------------

phylofile <- file.path("data/PhySpeTree","dist_bac_hcp.csv") #args[1]
chemfile <- file.path("data/processed","2018_Broc_ED09C.csv.mat.csv")#args[2]

PD_dist <- read.csv(file=phylofile)%>%
  mutate(ID=paste0(Sp1,"-",Sp2)) %>% 
  filter(PhyloDist!=0)

Chem_Dist <- read.csv(chemfile) %>%
  mutate(ID=paste0(Sp1,"-",Sp2)) %>% 
  filter(Chemdist!=0)

#Merge datasets.
df <- full_join(PD_dist,Chem_Dist,by="ID")  %>% 
  select(!c(Sp1.x,Sp2.x)) %>% 
  rename(Sp1=Sp1.y,Sp2=Sp2.y) %>% 
  filter(!is.na(Chemdist)) %>% 
  mutate(S1 = ifelse(Sp1 %in% c("ebw","ecr"),"E. coli",ifelse(Sp1 %in% c("pae","pau"),"P. aeruginosa",ifelse(Sp1 %in% c("seo","stm"),"S. enterica",ifelse(Sp1 =="cal","C. albicans",ifelse(Sp1 =="cne","C. neoformans",ifelse(Sp1=="sce","S. cerevisiae",ifelse(Sp1=="cgi","C. gattii",ifelse(Sp1=="spo","S. pombe","Unknown")))))))))%>% 
  mutate(S2 = ifelse(Sp2 %in% c("ebw","ecr"),"E. coli",ifelse(Sp2 %in% c("pae","pau"),"P. aeruginosa",ifelse(Sp2 %in% c("seo","stm"),"S. enterica",ifelse(Sp2 =="cal","C. albicans",ifelse(Sp2 =="cne","C. neoformans",ifelse(Sp2=="sce","S. cerevisiae",ifelse(Sp2=="cgi","C. gattii",ifelse(Sp2=="spo","S. pombe","Unknown"))))))))) %>% 
  mutate(class=ifelse(S1==S2,"Intraspecific comparison",ifelse(S1=="P. aeruginosa","Interspecific comparison (with Pseudomonas)",ifelse(S2=="P. aeruginosa","Interspecific comparison (with Pseudomonas)","Interspecific comparison (other than Pseudomonas)")))) %>% 
  mutate(class=fct_relevel(class,levels=c("Intraspecific comparison","Interspecific comparison (other than Pseudomonas)","Interspecific comparison (with Pseudomonas)")))

#Model statistics.
lmod <- lm(Chemdist~log(PhyloDist),data=df)
tid <- tidy(lmod)
aug <- augment(lmod)
gla <- glance(lmod)


g1 <- ggplot(df,aes(x=PhyloDist,y=Chemdist,group=class,label=ID))+
  geom_point(aes(shape=class))+
  geom_text_repel()+
  geom_smooth(method = "lm", se=FALSE, color="black")+
  labs(x="Phylogenetic Distance",
       y="Chemical Euclidean Distance",
       title="Chemical divergence accumulates with time",
       subtitle=paste("File:",chemfile),
       caption=paste("R^2=",round(gla$r.squared,4)))+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  annotation_logticks(sides="b")+
  theme_ipsum_rc(grid="XY")+
  geom_mark_circle(aes(label = class, fill=class,label.fontsize = 8))+
  theme_minimal()+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

d1 <- df %>% filter(class=="Intraspecific comparison")
#Model statistics.
lmod <- lm(Chemdist~log(PhyloDist),data=d1)
tid <- tidy(lmod)
aug <- augment(lmod)
gla <- glance(lmod)

d1 %>% 
  ggplot(aes(x=PhyloDist,y=Chemdist,group=class,label=ID))+
  geom_point(aes(shape=class))+
  geom_text_repel()+
  geom_smooth(method = "lm", se=FALSE, color="black")+
  labs(x="Phylogenetic Distance",
       y="Chemical Euclidean Distance",
       title="Chemical divergence accumulates with time",
       subtitle=paste("File:",chemfile),
       caption=paste("R^2=",round(gla$r.squared,4)))+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  annotation_logticks(sides="b")+
  theme_ipsum_rc(grid="XY")+
  theme_minimal()+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))


d1 <- df %>% filter(class=="Interspecific comparison (other than Pseudomonas)")
#Model statistics.
lmod <- lm(Chemdist~log(PhyloDist),data=d1)
tid <- tidy(lmod)
aug <- augment(lmod)
gla <- glance(lmod)

d1 %>% 
  ggplot(aes(x=PhyloDist,y=Chemdist,group=class,label=ID))+
  geom_point(aes(shape=class))+
  geom_text_repel()+
  geom_smooth(method = "lm", se=FALSE, color="black")+
  labs(x="Phylogenetic Distance",
       y="Chemical Euclidean Distance",
       title="Chemical divergence accumulates with time",
       subtitle=paste("File:",chemfile),
       caption=paste("R^2=",round(gla$r.squared,4)))+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  #annotation_logticks(sides="b")+
  theme_ipsum_rc(grid="XY")+
  theme_minimal()+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))


d1 <- df %>% filter(class=="Interspecific comparison (with Pseudomonas)")
#Model statistics.
lmod <- lm(Chemdist~log(PhyloDist),data=d1)
tid <- tidy(lmod)
aug <- augment(lmod)
gla <- glance(lmod)

d1 %>% 
  ggplot(aes(x=PhyloDist,y=Chemdist,group=class,label=ID))+
  geom_point(aes(shape=class))+
  geom_text_repel()+
  geom_smooth(method = "lm", se=FALSE, color="black")+
  labs(x="Phylogenetic Distance",
       y="Chemical Euclidean Distance",
       title="Chemical divergence accumulates with time",
       subtitle=paste("File:",chemfile),
       caption=paste("R^2=",round(gla$r.squared,4)))+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  #annotation_logticks(sides="b")+
  theme_ipsum_rc(grid="XY")+
  theme_minimal()+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))



ggplotly(g1)


ggplot(df,aes(x=PhyloDist,y=Chemdist,label=ID))+
  geom_point(aes(shape=class))+
  geom_text_repel()+
  geom_smooth(method = "lm", formula=y~log(x),se=TRUE, color="black")+
  labs(x="Phylogenetic Distance",
       y="Chemical Euclidean Distance",
       title="Chemical divergence accumulates with time",
       subtitle=paste("File:",chemfile))+
  theme_ipsum_rc(grid="XY")+
  geom_mark_circle(aes(label = class, fill=class,label.fontsize = 8))+
  theme_minimal()+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

ggsave(file.path("results",paste0("Fig_",args[2],".png")),plot=g1,device="png",width=15,height=10)
