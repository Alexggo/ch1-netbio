#' ---
#' title: "3.1.Heatmaps_PCA_morphospace"
#' author: "Alex Gil"
#' 
#' 
## -----------------------------------------------------------------------------
library(pacman)
p_load(tidyverse,RColorBrewer,plotly,dendextend,ggrepel,hrbrthemes,
       ggforce,ComplexHeatmap,ggbiplot,factoextra,OUwie,
       ouch,geiger,ape,phangorn,phytools,treeio,ggtree)




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
#' 
## -----------------------------------------------------------------------------

# Read data matrix
file <- file.path("data/4.PhylogeneticComparativeMethods/","dataset1_ratios.csv")

#Remove DDI with at least one NA.
dataset <- read_csv(file) %>%
        filter(!is.na(ebw)&!is.na(ecr)&!is.na(seo)&!is.na(stm)&!is.na(pae)&!is.na(pau))

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


mat_red <- dataset %>% select(ebw,ecr,stm,seo,pae,pau) %>%as.matrix
rownames(mat_red) <- dataset$drug_pair

dend <- hclust(dist(mat_red,method="euclidean"),method="average")  

annot <- dataset %>% select(drug_pair,drug_category1,drug_category2,categorycategory,
                            targeted_cellular_process1,targeted_cellular_process2,processprocess,
                            use1,use2,useuse,
                            drug_category,targeted_process,use,clusters,
                            sigma.rate)   %>% 
        as.matrix()
rownames(annot) <- annot[,1]


row_right_ha1 = rowAnnotation(category=annot[,11],
                              col = list(category = c("Same" = "red", "Different" = "black")))
row_right_ha2 = rowAnnotation(process=annot[,12],
                              col = list(process = c("Same" = "red", "Different" = "black")))
row_right_ha3 = rowAnnotation(use=annot[,13],
                              col = list(use = c("Same" = "red", "Different" = "black")))
row_right_ha4 = rowAnnotation(rate=anno_points((as.numeric(annot[,15]))))


row_left_ha = rowAnnotation(cluster=anno_block(gpar(fill = 1:13),
                            labels=c("13","5","3","11","7","6",
                                     "9","10","4","8","2","12","1"),
                            labels_gp = gpar(col = "white", fontsize = 10)))

Heatmap(mat_red,show_column_names = T,
        show_row_names = F,
        column_dend_height = unit(1,"cm"),
        show_column_dend = T,row_title = "Drug-Drug Interactions",
        column_title_side = "top",column_title = "Strains",column_names_side = "top",
        row_names_gp = gpar(cex=0.75),row_dend_width = unit(2,"cm"),
        name = "Interaction score",col=rev(brewer.pal(10,'RdBu')),
        left_annotation = row_left_ha,
        row_split = annot[,14])+row_right_ha1+row_right_ha2+row_right_ha3+row_right_ha4


## Plot PCA DDI~strain.
## -----------------------------------------------------------------------------
dataset1 <- dataset %>% as.data.frame()
res.pca <- dataset1[,10:15] 
colnames(res.pca) <- c("E. coli-EBW","E. coli-ECR",
                       "S. enterica-SEO","S. enterica-STM",
                       "P.aeruginosa-PAE","P.aeruginosa-PAU")

pca <- prcomp(res.pca, center = TRUE, scale = TRUE)

fviz_pca_biplot(pca, label ="var",
                repel = TRUE,
                habillage=dataset1$clusters,
                addEllipses=TRUE, ellipse.level=0.95) +
        theme_minimal() # Individuals color
fviz_eig(pca)# Eigenvalues


## Plot PCA Strain~DDI
## -----------------------------------------------------------------------------
dataset1 <- dataset %>% as.data.frame()
res.pca <- dataset1[,10:15] %>% t()
rownames(res.pca) <- c("E. coli-EBW","E. coli-ECR",
                       "S. enterica-SEO","S. enterica-STM",
                       "P.aeruginosa-PAE","P.aeruginosa-PAU")
pca <- prcomp(res.pca, center = TRUE, scale = TRUE)

#PCA
summary(pca)
pca.x=pca$x
pca.1.2=pca.x[,c(1,2)]


tree_nw$tip.label <- c("E. coli-EBW","E. coli-ECR",
                       "S. enterica-SEO","S. enterica-STM",
                       "P.aeruginosa-PAE","P.aeruginosa-PAU")


#Phylomorphospace plots.
phylomorphospace(tree_nw, pca$x[,1:2], 
                 label = "horizontal", node.size=c(.5,1),
                 xlim=c(-50,60),
                 xlab="PC1 (47.01%)",
                 ylab="PC2 (21.54%)")


#phylomorphospace3d(tree_nw,pca.x[,c(1,2,3)],method="static")
fancyTree(tree_nw,X=pca.x[,c(1,2)],
          type="traitgram3d",method="static")


###############

## Same for filtered DDI

#Only same resistence
same_res <- c("Imipenem","Aztreonam","Cephalexin","Teicoplanin",
              "Bacitracin","Metronidazole","Clofazimine",
              "Doxorubicin","Mitomycin C","Ciprofloxacin","Levofloxacin",
              "Sulfamonomethoxine","Cerulenin","Nisin","Daptomycin",
              "Nonactin","Paraquat","Linezolid","Tigecycline",
              "Spiramycin","Erythromycin","Clarithromycin",
              "Caffeine","Reserpine","Procaine","Theophylline",
              "Verapamil","Meropenem","Triclosan",
              "Gentamicin","Berberine","Metformin","Diclofenac")

dataset_filt <- dataset |> 
  filter(Drug1 %in% same_res) |> 
  filter(Drug2 %in% same_res) |> 
  filter(!is.na(ebw)&!is.na(ecr)&!is.na(seo)&!is.na(stm)&!is.na(pae)&!is.na(pau))


dataset_filt |> dim()

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


mat_red <- dataset_filt %>% select(ebw,ecr,stm,seo,pae,pau) %>%as.matrix
rownames(mat_red) <- dataset$drug_pair

dend <- hclust(dist(mat_red,method="euclidean"),method="average")  

annot <- dataset %>% select(drug_pair,drug_category1,drug_category2,categorycategory,
                            targeted_cellular_process1,targeted_cellular_process2,processprocess,
                            use1,use2,useuse,
                            drug_category,targeted_process,use,clusters,
                            sigma.rate)   %>% 
  as.matrix()
rownames(annot) <- annot[,1]


row_right_ha1 = rowAnnotation(category=annot[,11],
                              col = list(category = c("Same" = "red", "Different" = "black")))
row_right_ha2 = rowAnnotation(process=annot[,12],
                              col = list(process = c("Same" = "red", "Different" = "black")))
row_right_ha3 = rowAnnotation(use=annot[,13],
                              col = list(use = c("Same" = "red", "Different" = "black")))
row_right_ha4 = rowAnnotation(rate=anno_points((as.numeric(annot[,15]))))


row_left_ha = rowAnnotation(cluster=anno_block(gpar(fill = 1:13),
                                               labels=c("13","5","3","11","7","6",
                                                        "9","10","4","8","2","12","1"),
                                               labels_gp = gpar(col = "white", fontsize = 10)))

Heatmap(mat_red,show_column_names = T,
        show_row_names = F,
        column_dend_height = unit(1,"cm"),
        show_column_dend = T,row_title = "Drug-Drug Interactions",
        column_title_side = "top",column_title = "Strains",column_names_side = "top",
        row_names_gp = gpar(cex=0.75),row_dend_width = unit(2,"cm"),
        name = "Interaction score",col=rev(brewer.pal(10,'RdBu')),
        left_annotation = row_left_ha,
        row_split = annot[,14])+row_right_ha1+row_right_ha2+row_right_ha3+row_right_ha4


## Plot PCA DDI~strain.
## -----------------------------------------------------------------------------
dataset1 <- dataset %>% as.data.frame()
res.pca <- dataset1[,10:15] 
colnames(res.pca) <- c("E. coli-EBW","E. coli-ECR",
                       "S. enterica-SEO","S. enterica-STM",
                       "P.aeruginosa-PAE","P.aeruginosa-PAU")

pca <- prcomp(res.pca, center = TRUE, scale = TRUE)

fviz_pca_biplot(pca, label ="var",
                repel = TRUE,
                habillage=dataset1$clusters,
                addEllipses=TRUE, ellipse.level=0.95) +
  theme_minimal() # Individuals color
fviz_eig(pca)# Eigenvalues


## Plot PCA Strain~DDI
## -----------------------------------------------------------------------------
dataset1 <- dataset %>% as.data.frame()
res.pca <- dataset1[,10:15] %>% t()
rownames(res.pca) <- c("E. coli-EBW","E. coli-ECR",
                       "S. enterica-SEO","S. enterica-STM",
                       "P.aeruginosa-PAE","P.aeruginosa-PAU")
pca <- prcomp(res.pca, center = TRUE, scale = TRUE)

#PCA
summary(pca)
pca.x=pca$x
pca.1.2=pca.x[,c(1,2)]


tree_nw$tip.label <- c("E. coli-EBW","E. coli-ECR",
                       "S. enterica-SEO","S. enterica-STM",
                       "P.aeruginosa-PAE","P.aeruginosa-PAU")


#Phylomorphospace plots.
phylomorphospace(tree_nw, pca$x[,1:2], 
                 label = "horizontal", node.size=c(.5,1),
                 xlim=c(-50,60),
                 xlab="PC1 (47.01%)",
                 ylab="PC2 (21.54%)")


#phylomorphospace3d(tree_nw,pca.x[,c(1,2,3)],method="static")
fancyTree(tree_nw,X=pca.x[,c(1,2)],
          type="traitgram3d",method="static")




## Same but only same Res/Sen


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
#' 
## -----------------------------------------------------------------------------

# Read data matrix
file <- file.path("data/B.SameRes_across_strains_DDIs/4.PhylogeneticComparativeMethods/","dataset1_ratios.csv")

#Remove DDI with at least one NA.
dataset <- read_csv(file) %>%
  filter(!is.na(ebw)&!is.na(ecr)&!is.na(seo)&!is.na(stm)&!is.na(pae)&!is.na(pau))

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


mat_red <- dataset %>% select(ebw,ecr,stm,seo,pae,pau) %>%as.matrix
rownames(mat_red) <- dataset$drug_pair

dend <- hclust(dist(mat_red,method="euclidean"),method="average")  

annot <- dataset %>% select(drug_pair,drug_category1,drug_category2,categorycategory,
                            targeted_cellular_process1,targeted_cellular_process2,processprocess,
                            use1,use2,useuse,
                            drug_category,targeted_process,use,clusters,
                            sigma.rate)   %>% 
  as.matrix()
rownames(annot) <- annot[,1]


row_right_ha1 = rowAnnotation(category=annot[,11],
                              col = list(category = c("Same" = "red", "Different" = "black")))
row_right_ha2 = rowAnnotation(process=annot[,12],
                              col = list(process = c("Same" = "red", "Different" = "black")))
row_right_ha3 = rowAnnotation(use=annot[,13],
                              col = list(use = c("Same" = "red", "Different" = "black")))
row_right_ha4 = rowAnnotation(rate=anno_points((as.numeric(annot[,15]))))


row_left_ha = rowAnnotation(cluster=anno_block(gpar(fill = 1:13),
                                               labels=c("13","5","3","11","7","6",
                                                        "9","10","4","8","2","12","1"),
                                               labels_gp = gpar(col = "white", fontsize = 10)))

Heatmap(mat_red,show_column_names = T,
        show_row_names = F,
        column_dend_height = unit(1,"cm"),
        show_column_dend = T,row_title = "Drug-Drug Interactions",
        column_title_side = "top",column_title = "Strains",column_names_side = "top",
        row_names_gp = gpar(cex=0.75),row_dend_width = unit(2,"cm"),
        name = "Interaction score",col=rev(brewer.pal(10,'RdBu')),
        left_annotation = row_left_ha,
        row_split = annot[,14])+row_right_ha1+row_right_ha2+row_right_ha3+row_right_ha4


## Plot PCA DDI~strain.
## -----------------------------------------------------------------------------
dataset1 <- dataset %>% as.data.frame()
res.pca <- dataset1[,12:17] 
colnames(res.pca) <- c("E. coli-EBW","E. coli-ECR",
                       "S. enterica-SEO","S. enterica-STM",
                       "P.aeruginosa-PAE","P.aeruginosa-PAU")

pca <- prcomp(res.pca, center = TRUE, scale = TRUE)

fviz_pca_biplot(pca, label ="var",
                repel = TRUE,
                habillage=dataset1$clusters,
                addEllipses=TRUE, ellipse.level=0.95) +
  theme_minimal() # Individuals color
fviz_eig(pca)# Eigenvalues


## Plot PCA Strain~DDI
## -----------------------------------------------------------------------------
dataset1 <- dataset %>% as.data.frame()
res.pca <- dataset1[,12:17] %>% t()
rownames(res.pca) <- c("E. coli-EBW","E. coli-ECR",
                       "S. enterica-SEO","S. enterica-STM",
                       "P.aeruginosa-PAE","P.aeruginosa-PAU")
pca <- prcomp(res.pca, center = TRUE, scale = TRUE)

#PCA
summary(pca)
pca.x=pca$x
pca.1.2=pca.x[,c(1,2)]


tree_nw$tip.label <- c("E. coli-EBW","E. coli-ECR",
                       "S. enterica-SEO","S. enterica-STM",
                       "P.aeruginosa-PAE","P.aeruginosa-PAU")


#Phylomorphospace plots.
phylomorphospace(tree_nw, pca$x[,1:2], 
                 label = "horizontal", node.size=c(.5,1),
                 xlim=c(-50,60),
                 xlab="PC1 (47.01%)",
                 ylab="PC2 (21.54%)")


#phylomorphospace3d(tree_nw,pca.x[,c(1,2,3)],method="static")
fancyTree(tree_nw,X=pca.x[,c(1,2)],
          type="traitgram3d",method="static")
