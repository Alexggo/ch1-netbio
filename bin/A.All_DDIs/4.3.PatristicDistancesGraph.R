#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#' ---
#' title: 4.3.PatristicDistancesGraph.
#' author: "Alex Gil"
#' 
#' 
## -----------------------------------------------------------------------------
library(pacman)
p_load(plyr,tidyverse,broom,ape,phangorn,phytools,treeio,ggtree,plotly,
dendextend,RColorBrewer,car,factoextra,OUwie,ouch,factoextra,geiger,ggrepel,
hrbrthemes,ggforce,ComplexHeatmap,ggbiplot,patchwork,performance)



#' Input files, phylogenetic tree and matrix.
#' 
#'   -Choose the tree: args[1]
#'     + bac_species.txt.hcp
#'     + bac_species.txt.rna
#'     + fun_species.txt.hcp
#'     + fun_species.txt.rna
#'   -Choose the matrix (rows are species, columns are traits): args[2]
#'     + 2011_Spit_S1.1.csv
#'     + 2011_Spit_S3.1.csv
#'     + 2015_Spit_S2.1.csv
#'     + 2015_Spit_S5.1.csv
#'     + 2018_Broc_ED02.1.csv
#'     + 2018_Broc_ED03.2.csv
#'     + 2018_Broc_ED09C.csv.mat.csv
#' 
## -----------------------------------------------------------------------------

mat <- read.csv("data/1.processed/Broc2018_maindataset.csv") %>%
  filter(!is.na(ebw)&!is.na(ecr)&!is.na(seo)&!is.na(stm)&!is.na(pae)&!is.na(pau))

clust <- read.csv("data/3.InteractionScores_tSNE/tSNE_Clustermembership_ppx706.csv")
mat1 <- full_join(mat,clust,by="drug_pair")


# Read phylogenetic tree
# Read phylogenetic tree
species <- c("Escherichia_coli_K-12_ebw",
             "Escherichia_coli_O8_IAI1_ecr",
             "Salmonella_enterica_serovar_Typhimurium_LT2_stm",
             "Salmonella_enterica_serovar_Typhimurium_14028S_seo",
             "Pseudomonas_aeruginosa_PAO1_pae",
             "Pseudomonas_aeruginosa_UCBPP-PA14_pau")

treetype <- "bac_species.txt.hcp"
treefile <- file.path('data/2.2.Phylogenetics_Bayesian/DDI_BD_str',"concatenate.trees.tre")
tree_nw <- read.nexus(treefile)
# Drop the species that are not needed.
not.species <- tree_nw$tip.label[!(tree_nw$tip.label %in% species)]
tree_nw <- drop.tip(tree_nw,not.species)
tree_nw$tip.label <- c("ebw","ecr","pae","pau","seo","stm")

PatristicDistMatrix<-cophenetic(tree_nw) |> 
  as.data.frame()
PatristicDistMatrix$Sp1 <- rownames(PatristicDistMatrix)

PD_dist <- PatristicDistMatrix |> pivot_longer(names_to = "Sp2",
                                    values_to="PhyloDist",
                                    1:6) |>
  mutate(ID=paste0(Sp1,"-",Sp2)) %>% 
  filter(PhyloDist!=0) |> 
  mutate(test=Sp1<Sp2) |> 
  filter(test==TRUE) |> 
  select(-test)

#Euclidean chemical distances.
g <- mat1 %>% select(12:17) %>% t()
dend <- hclust(dist(g,method="euclidean"),method="average")  %>% 
  as.phylo()
PatristicDistMatrix<-cophenetic(dend)
tabnew<- PatristicDistMatrix %>% as.data.frame()
tabnew <- tabnew[ , order(names(tabnew))]
tabnew$Sp1 <- rownames(tabnew)
tabnew <- tabnew %>%
    arrange(Sp1) %>%
    select(-Sp1)

tabnew[lower.tri(tabnew)] <- NA
tabnew$Sp1 <- colnames(tabnew)
tabnew <- tabnew %>%
    select(dim(tabnew)[2],1:(dim(tabnew)[2]-1))
tabnew <- tabnew  %>%
    pivot_longer(names_to="Sp2",values_to="Chemdist",2:(dim(tabnew)[2]))%>%
    filter(!is.na(Chemdist)) %>%
    filter(Chemdist!=0) %>%
    mutate(ID=paste0(Sp1,"-",Sp2))
Chem_Dist <- tabnew

df <- inner_join(PD_dist,Chem_Dist,by="ID")  %>% 
  mutate(Sp1=Sp1.y,Sp2=Sp2.y) %>% 
  select(!c(Sp1.x,Sp2.x,Sp1.y,Sp2.y)) %>% 
  filter(!is.na(Chemdist)) %>% 
  mutate(S1 = ifelse(Sp1 %in% c("ebw","ecr"),"E. coli",ifelse(Sp1 %in% c("pae","pau"),"P. aeruginosa",ifelse(Sp1 %in% c("seo","stm"),"S. enterica",ifelse(Sp1 =="cal","C. albicans",ifelse(Sp1 =="cne","C. neoformans",ifelse(Sp1=="sce","S. cerevisiae",ifelse(Sp1=="cgi","C. gattii",ifelse(Sp1=="spo","S. pombe","Unknown")))))))))%>% 
  mutate(S2 = ifelse(Sp2 %in% c("ebw","ecr"),"E. coli",ifelse(Sp2 %in% c("pae","pau"),"P. aeruginosa",ifelse(Sp2 %in% c("seo","stm"),"S. enterica",ifelse(Sp2 =="cal","C. albicans",ifelse(Sp2 =="cne","C. neoformans",ifelse(Sp2=="sce","S. cerevisiae",ifelse(Sp2=="cgi","C. gattii",ifelse(Sp2=="spo","S. pombe","Unknown"))))))))) %>% 
  mutate(class=ifelse(S1==S2,"Intraspecific comparison",ifelse(S1=="P. aeruginosa","Interspecific comparison (with Pseudomonas)",ifelse(S2=="P. aeruginosa","Interspecific comparison (with Pseudomonas)","Interspecific comparison (other than Pseudomonas)"))))

#Model statistics.
lmod1 <- lm(Chemdist~PhyloDist,data=df)
gla1 <- glance(lmod1)
check_model(lmod1)

lmod2 <- lm(Chemdist~log(PhyloDist),data=df)
tid2 <- tidy(lmod2)
aug2 <- augment(lmod2)
gla2 <- glance(lmod2)
performance::check_model(lmod2)
check_outliers(lmod2)

g1 <- ggplot(df,aes(x=PhyloDist,y=Chemdist,label=ID))+
  geom_point(aes(color=class))+
  geom_text_repel()+
  geom_smooth(method = "lm", formula=y~log(x),se=FALSE, color="black")+
  labs(x="Phylogenetic Distance MRCA (mya)",
       y="Chemical Patristic Distance from Dendrogram")+
  theme_ipsum_rc(grid="XY")+
  theme_minimal()+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(legend.position = "none")+
  geom_mark_circle(aes(fill = class,label=class))

g2 <- ggplot(df,aes(x=PhyloDist,y=Chemdist,label=ID))+
  geom_point(aes(color=class))+
  geom_text_repel()+
  geom_smooth(method = "lm", se=FALSE, color="black")+
  labs(x="Phylogenetic Distance MRCA (mya)",
       y="Chemical Patristic Distance from Dendrogram")+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  annotation_logticks(sides="b")+
  theme_ipsum_rc(grid="XY")+
  theme_minimal()+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(legend.position = "none")+
  geom_mark_circle(aes(fill = class,label=class))+
  annotate("text", x=500, y=5.5, label= "Chemdist ~ log(PhyloDist)")+
  annotate("text", x=500, y=5.4, label= paste0("R squared:",round(gla2$r.squared,2),", p-value:",round(gla2$p.value,10)))


g1+g2
