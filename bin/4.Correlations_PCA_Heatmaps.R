#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#' ---
#' title: 4.1.PatristicDistancesGraph.
#' author: "Alex Gil"
#' 
#' 
## -----------------------------------------------------------------------------
library(pacman)
p_load(tidyverse, broom, ape, phangorn, phytools, treeio, ggtree, 
       plotly, dendextend, RColorBrewer, car, factoextra, OUwie, 
       ouch, geiger, ggrepel, hrbrthemes, ggforce, ComplexHeatmap, 
       ggbiplot, patchwork, performance,pdftools,png,ggplotify)


set_name <- "allddi" #allddi sen2in1
full_df <- read.csv(paste0("results/",set_name,"/","DDI_table_rates_",set_name,".csv"))


# Read phylogenetic tree
species <- c("Escherichia_coli_K-12_ebw",
             "Escherichia_coli_O8_IAI1_ecr",
             "Salmonella_enterica_serovar_Typhimurium_LT2_stm",
             "Salmonella_enterica_serovar_Typhimurium_14028S_seo",
             "Pseudomonas_aeruginosa_PAO1_pae",
             "Pseudomonas_aeruginosa_UCBPP-PA14_pau")

# Calculate distances from single tree.
treefile <- file.path('data/2.2.Phylogenetics_Bayesian',"Phylo-p1.nwk.tree")
tree_nw <- read.nexus(treefile)

tree1 <- tree_nw |> 
  ape::rotateConstr(c("stm","seo","ebw","ecr","pae","pau")) |> 
  ggtree(ladderize = FALSE,size=1.25)+geom_tiplab(size=7) +scale_x_continuous(labels=abs)

fig1.a <- revts(tree1)+
  theme_tree2()+ 
  xlab("Phylogenetic Distance (MYA)")+
  theme(axis.title = element_text(size = 16),  # Set axis title size
        axis.text = element_text(size = 15))

fig1.a
ggsave(paste0("results/", set_name, "/fig1A_sptree_", set_name, ".pdf"), 
       plot = fig1.a,bg='transparent',
       height = 15,
       width = 6)


### HEATMAP
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
  mutate(Category=ifelse(drug_category1==drug_category2,"Same","Different"),
         Use=ifelse(use1==use2,"Same","Different"),
         Process=ifelse(targeted_cellular_process1==targeted_cellular_process2,"Same","Different")) |> 
  as.matrix() |> 
  t()
colnames(annot) <- annot[1,]
annot[11,] <- as.character(as.numeric(annot[11,]))
annot[12,] <- as.character(as.numeric(annot[12,]))

anno_df = data.frame(
  Category=annot[15,],
  Process=annot[17,],
  Use=annot[16,]
)

Rates <- as.numeric(annot[12,])

ha = HeatmapAnnotation(df = anno_df,
                       col = list(Rates,
                                  Category = c("Same" = "green", "Different" = "black"),
                                  Process = c("Same" = "green", "Different" = "black"),
                                  Use = c("Same" = "green", "Different" = "black"))
                       
)


bottom= HeatmapAnnotation(`Evolutionary rate`=anno_barplot(Rates))

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


# Check colors with labels
Heatmap(mat_red,
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
        top_annotation = HeatmapAnnotation(df = data.frame(annot[11,])),
        bottom_annotation = h)

fig1.b = grid.grabExpr(draw(h1)) 

fig1.ab <- wrap_plots(fig1.a,fig1.b,widths = c(0.2,0.8))

a <- 3.5
pdf(paste0("results/",set_name,"/","fig1B_heatmap_",set_name,".pdf"),
    width = 3*a,height = 2*a)
h1
dev.off()

# Phylomorphospace and PCA

dataset1 <- full_df %>% as.data.frame()
res.pca <- dataset1[,19:24] 
colnames(res.pca) <- c("E. coli-EBW","E. coli-ECR",
                       "S. enterica-SEO","S. enterica-STM",
                       "P.aeruginosa-PAE","P.aeruginosa-PAU")

pca <- prcomp(res.pca, center = TRUE, scale = TRUE)

pca_fig1 <- fviz_pca_biplot(pca, label ="var",
                            repel = TRUE,
                            habillage=dataset1$clusters,
                            addEllipses=TRUE, ellipse.level=0.95) +
  theme_minimal() # Individuals color
pca_fig2 <- fviz_eig(pca)# Eigenvalues

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
phylomorphospace3d(tree_nw,pca.x[,c(1,2,3)],method="static")
fancyTree(tree_nw,X=pca.x[,c(1,2)],
          type="traitgram3d",method="static")
pca_fig3 <- recordPlot()
plot.new() ## clean up device

#This one is default:
phylomorphospace(tree_nw, pca$x[,1:2], 
                 label = "horizontal", node.size=c(.5,1),
                 xlim=c(-50,60),
                 xlab="PC1 (47.01%)",
                 ylab="PC2 (21.54%)")

obj <- phylomorphospace(tree_nw, pca$x[,1:2], 
                        label = "horizontal", node.size=c(.5,1),
                        xlim=c(-50,60),
                        xlab="PC1 (47.01%)",
                        ylab="PC2 (21.54%)")
grid()

Data<-data.frame(
  xstart=obj$xx[obj$edge[,1]],
  ystart=obj$yy[obj$edge[,1]],
  xstop=obj$xx[obj$edge[,2]],
  ystop=obj$yy[obj$edge[,2]],
  nodestart=obj$edge[,1],
  nodestop=obj$edge[,2],
  strains=c("","","ebw","ecr","","seo","stm","","pae","pau"))

fig1.c <- ggplot()+
  geom_point(data=Data[1,],aes(x = xstart, y = ystart), color = "red", size = 3,
             shape=8) +  # Points where segments start
  geom_segment(data=Data,aes(x=xstart,y=ystart,xend=xstop,yend=ystop), 
               size = 0.5,color="black")+
  geom_point(data=Data,aes(x = xstop, y = ystop), color = "black", size = 3)+     # Points where segments end
  theme_bw()+
  xlab('PC1 (47.01%)')+
  ylab('PC2 (21.54%)')+
  geom_text(data=Data,aes(x = xstop, y = ystop, label = strains), vjust = -1)+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        axis.title = element_text(size = 12),  # Set axis title size
        axis.text = element_text(size = 11))

ggsave(
  filename = paste0("results/allddi/fig1C_phylomorph.pdf"),
  plot = fig1.c,
  width = 7,
  height = 7
)

# Save figures
ggsave(
  filename = paste0("results/allddi/pca_fig1.pdf"),
  plot = pca_fig1,
  width = 15,
  height = 15
)

ggsave(
  filename = paste0("results/allddi/pca_fig2.pdf"),
  plot = pca_fig2,
  width = 7,
  height = 7
)

pdf(paste0("results/allddi/pca_fig3.pdf"))
pca_fig3 # redraw
dev.off()

# Correlations
PatristicDistMatrix<-cophenetic(tree_nw)
PatristicDistMatrix <- PatristicDistMatrix/2

PatristicDistMatrix <- PatristicDistMatrix|> 
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
g <- full_df %>% select(ebw,ecr,seo,stm,pae,pau) %>% t()
dend <- hclust(dist(g,method="euclidean"),method="average")  %>% 
  as.phylo()

dend2 <- rotate(dend,7)
dend3 <- rotate(dend2,9)
dend4 <- rotate(dend3,11)

plot(dend4)
add.scale.bar()

PatristicDistMatrix<-cophenetic(dend)
PatristicDistMatrix <- PatristicDistMatrix/2
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

fig1.d<- ggplot(df,aes(x=PhyloDist,y=Chemdist,label=ID))+
  geom_point(aes(color=class))+
  geom_text_repel()+
  geom_smooth(method = "lm", formula=y~log(x),se=FALSE, color="black")+
  labs(x="Phylogenetic Distance (MYA)",
       y="Chemical Patristic Distance from Dendrogram")+
  theme_ipsum_rc(grid="XY")+
  theme_minimal()+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        axis.title = element_text(size = 15),  # Set axis title size
        axis.text = element_text(size = 13))+
  geom_mark_circle(aes(fill = class,label=class))+
  coord_cartesian(xlim = c(0, 1550))  # Set x-axis limits


ggsave(
  filename = paste0("results/", set_name, "/fig1D_dist_", set_name, ".pdf"),
  plot = fig1.d,
  height = 7,
  width = 7
)

lmod2 <- lm(Chemdist~log(PhyloDist),data=df)
tid2 <- tidy(lmod2)
aug2 <- augment(lmod2)
gla2 <- glance(lmod2)

supfig.1a <- ggplot(df,aes(x=PhyloDist,y=Chemdist,label=ID))+
  geom_point(aes(color=class))+
  geom_text_repel()+
  geom_smooth(method = "lm", se=FALSE, color="black")+
  labs(x="Phylogenetic Distance (MYA)",
       y="Chemical Patristic Distance from Dendrogram")+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  annotation_logticks(sides="b")+
  theme_ipsum_rc(grid="XY")+
  theme_minimal()+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        axis.title = element_text(size = 15),  # Set axis title size
        axis.text = element_text(size = 13))+
  geom_mark_circle(aes(fill = class,label=class))+
  annotate("text", x=60, y=2.65, label= "ChemDist ~ log(PhyloDist)",
           size=4.5,hjust=0)+
  annotate("text", x=60, y=2.58, label= paste0("R squared:",round(gla2$r.squared,2),", p-value:",round(gla2$p.value,10)),
           size=4.5,hjust=0)

ggsave(
  filename = paste0("results/", set_name, "/supfig_logcorr_", set_name, ".pdf"),
  plot = supfig.1a,
  height = 7,
  width = 7
)

ggsave(
  filename = paste0("results/", set_name, "/supfig_logcorr_", set_name, ".png"),
  plot = supfig.1a,
  height = 7,
  width = 7,dpi=600
)


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

mat <- full_df |> 
  filter(Drug1 %in% same_res) |> 
  filter(Drug2 %in% same_res) |> 
  filter(!is.na(ebw)&!is.na(ecr)&!is.na(seo)&!is.na(stm)&!is.na(pae)&!is.na(pau))


PatristicDistMatrix<-cophenetic(tree_nw)
PatristicDistMatrix <- PatristicDistMatrix/2

PatristicDistMatrix <- PatristicDistMatrix|> 
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
g <- mat %>% select(ebw,ecr,seo,stm,pae,pau) %>% t()
dend <- hclust(dist(g,method="euclidean"),method="average")  %>% 
  as.phylo()

PatristicDistMatrix<-cophenetic(dend)
PatristicDistMatrix <- PatristicDistMatrix/2

tabnew <- PatristicDistMatrix|> 
  as.data.frame()

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

lmod2 <- lm(Chemdist~log(PhyloDist),data=df)
tid2 <- tidy(lmod2)
aug2 <- augment(lmod2)
gla2 <- glance(lmod2)

g1 <- ggplot(df,aes(x=PhyloDist,y=Chemdist,label=ID))+
  geom_point(aes(color=class))+
  geom_text_repel()+
  geom_smooth(method = "lm", formula=y~log(x),se=FALSE, color="black")+
  labs(x="Phylogenetic Distance (MYA)",
       y="Chemical Patristic Distance from Dendrogram")+
  theme_ipsum_rc(grid="XY")+
  theme_minimal()+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA))+
  theme(legend.position = "none")+
  geom_mark_circle(aes(fill = class,label=class))

g2 <- ggplot(df,aes(x=PhyloDist,y=Chemdist,label=ID))+
  geom_point(aes(color=class))+
  geom_text_repel()+
  geom_smooth(method = "lm", se=FALSE, color="black")+
  labs(x="Phylogenetic Distance (MYA)",
       y="Chemical Patristic Distance from Dendrogram")+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  annotation_logticks(sides="b")+
  theme_ipsum_rc(grid="XY")+
  theme_minimal()+
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA))+
  theme(legend.position = "none")+
  geom_mark_circle(aes(fill = class,label=class))+
  annotate("text", x=400, y=0.82, label= "Chemdist ~ log(PhyloDist)")+
  annotate("text", x=400, y=0.8, label= paste0("R squared:",round(gla2$r.squared,2),", p-value:",round(gla2$p.value,10)))

ggsave(
  filename = paste0("results/onlysen", "/fig_dist_onlysen.pdf"),
  plot = g1,
  height = 7,
  width = 7
)

ggsave(
  filename = paste0("results/onlysen", "/fig_dist_onlysen.pdf"),
  plot = g2,
  height = 7,
  width = 7
)

