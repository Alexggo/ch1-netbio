#' ---
#' title: "4.2.Heatmaps_PCA_morphospace"
#' author: "Alex Gil"
#' 
#' 
## -----------------------------------------------------------------------------
library(pacman)
p_load(tidyverse,RColorBrewer,plotly,dendextend,ggrepel,hrbrthemes,
       ggforce,ggbiplot,factoextra,OUwie,
       ouch,geiger,ape,phangorn,phytools,treeio,ggtree,ComplexHeatmap)

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

set_name <- "sen2in1" #sen2in1 allddi
full_df <- read.csv(paste0("results/",set_name,"/","DDI_table_rates_",set_name,".csv"))

###
# Read phylogenetic tree
species <- c("Escherichia_coli_K-12_ebw",
             "Escherichia_coli_O8_IAI1_ecr",
             "Salmonella_enterica_serovar_Typhimurium_LT2_stm",
             "Salmonella_enterica_serovar_Typhimurium_14028S_seo",
             "Pseudomonas_aeruginosa_PAO1_pae",
             "Pseudomonas_aeruginosa_UCBPP-PA14_pau")

treefile <- file.path('data/2.2.Phylogenetics_Bayesian',"Phylo-p1.nwk.tree")
tree_nw <- read.nexus(treefile)
# Drop the species that are not needed.
not.species <- tree_nw$tip.label[!(tree_nw$tip.label %in% species)]
tree_nw <- drop.tip(tree_nw,not.species)

tree_nw$tip.label <- c("ebw","ecr","pae","pau","seo","stm")
####
sum_df <- full_df  |> 
  group_by(clusters,sigma.rate) |> 
  dplyr::summarize(n = n())|> 
  arrange(desc(sigma.rate))

mat_red <- full_df %>% select(ebw,ecr,stm,seo,pae,pau) %>%as.matrix |> t()
colnames(mat_red) <- full_df$drug_pair

annot <- full_df |> 
  select(drug_pair,drug_category1,drug_category2,categorycategory,
                            targeted_cellular_process1,targeted_cellular_process2,processprocess,
                            use1,use2,useuse,
                            clusters,
                            sigma.rate,
                            SR_score_total,sum_g)    |>  
  mutate(category=ifelse(drug_category1==drug_category2,"Same","Different"),
         use=ifelse(use1==use2,"Same","Different"),
         process=ifelse(targeted_cellular_process1==targeted_cellular_process2,"Same","Different")) |> 
  as.matrix() |> 
  t()
colnames(annot) <- annot[1,]
annot[11,] <- as.character(as.numeric(annot[11,]))
annot[12,] <- as.character(as.numeric(annot[12,]))

anno_df = data.frame(
  category=annot[15,],
  process=annot[17,],
  use=annot[16,]
)

rates <- as.numeric(annot[12,])

ha = HeatmapAnnotation(df = anno_df,
                       col = list(rates,
                                  category = c("Same" = "green", "Different" = "black"),
                                  process = c("Same" = "green", "Different" = "black"),
                                  use = c("Same" = "green", "Different" = "black"))
                       
)


bottom= HeatmapAnnotation(rate=anno_barplot(rates))

h=c(ha,bottom)

cl_1 <- brewer.pal(9,"Reds")[1:(length(sum_df$clusters)/2)]
cl_2 <- brewer.pal(9,"Greens")[1:(length(sum_df$clusters)/2)]
cl_v <- c(cl_1,cl_2)
cl_v <- cl_v[1:length(sum_df$clusters)]
names(cl_v) <- sum_df$clusters
t1 = HeatmapAnnotation(df = data.frame(clusters=annot[11,]),
                      col=list(clusters=cl_v))


list_clusters <- list(allddi=as.character(c(14,7,3,4,16,8,5,13,10,15,9,12,6,2,11,1)),
                      sen2in1 = as.character(c(7,11,9,10,3,8,4,6,5,2,1)))
t2 <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = 1:length(sum_df$clusters)),
                                             labels = list_clusters[[set_name]],
                                             labels_gp = gpar(col = "white", fontsize = 14)))

h1 <- Heatmap(mat_red,
        show_column_names = F,
        show_row_names = T,row_names_side = "left",
        column_dend_height = unit(3,"cm"),
        show_column_dend = T,
        row_title = "Strains",
        column_title_side = "top",
        column_title = "Drug-Drug Interactions",
        column_names_side = "top",
        row_names_gp = gpar(cex=0.75),
        row_dend_width = unit(1,"cm"),
        name = "Interaction score",
        col=rev(brewer.pal(10,'RdBu')),
        column_split = annot[11,],
        show_heatmap_legend = T,
        top_annotation = t2,
        bottom_annotation = h)

## Plot PCA DDI~strain.
## -----------------------------------------------------------------------------
dataset1 <- full_df %>% as.data.frame()
res.pca <- dataset1[,19:24] 
colnames(res.pca) <- c("E. coli-EBW","E. coli-ECR",
                       "S. enterica-SEO","S. enterica-STM",
                       "P.aeruginosa-PAE","P.aeruginosa-PAU")

pca <- prcomp(res.pca, center = TRUE, scale = TRUE)

p1 <- fviz_pca_biplot(pca, label ="var",
                repel = TRUE,
                habillage=dataset1$clusters,
                addEllipses=TRUE, ellipse.level=0.95) +
        theme_minimal() # Individuals color
p2 <- fviz_eig(pca)# Eigenvalues


## Plot PCA Strain~DDI
## -----------------------------------------------------------------------------
dataset1 <- full_df %>% as.data.frame()
res.pca <- dataset1[,19:24] %>% t()
rownames(res.pca) <- c("ebw","ecr",
                       "seo","stm",
                       "pae","pau")
pca <- prcomp(res.pca, center = TRUE, scale = TRUE)

#PCA
summary(pca)
pca.x=pca$x
pca.1.2=pca.x[,c(1,2)]


#Phylomorphospace plots.
phylomorphospace(tree_nw, pca$x[,1:2], 
                 label = "horizontal", node.size=c(.5,1),
                 xlim=c(-50,60),
                 xlab="PC1 (47.01%)",
                 ylab="PC2 (21.54%)")


#phylomorphospace3d(tree_nw,pca.x[,c(1,2,3)],method="static")
fancyTree(tree_nw,X=pca.x[,c(1,2)],
          type="traitgram3d",method="static")


a <- 3
pdf(paste0("results/",set_name,"/","fig1B_heatmap_",set_name,".pdf"),
    width = 3*a,height = 2*a)
h1
dev.off()

pdf(paste0("results/",set_name,"/","fig1C_pca_",set_name,".pdf"),
    height=5,width = 5)
p1
p2

#Phylomorphospace plots.
phylomorphospace(tree_nw, pca$x[,1:2], 
                 label = "horizontal", node.size=c(.5,1),
                 xlim=c(-50,60),
                 xlab="PC1 (47.01%)",
                 ylab="PC2 (21.54%)")


#phylomorphospace3d(tree_nw,pca.x[,c(1,2,3)],method="static")
fancyTree(tree_nw,X=pca.x[,c(1,2)],
          type="traitgram3d",method="static")

dev.off()
