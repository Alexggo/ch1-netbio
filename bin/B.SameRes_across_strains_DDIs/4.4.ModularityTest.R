library(pacman)
p_load(geomorph,tidyverse,broom,ape)


mat <- read.csv("data/1.processed/Broc2018_maindataset.csv") %>%
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


same_res <- c("Imipenem","Aztreonam","Cephalexin","Teicoplanin",
              "Bacitracin","Metronidazole","Clofazimine",
              "Doxorubicin","Mitomycin C","Ciprofloxacin","Levofloxacin",
              "Sulfamonomethoxine","Cerulenin","Nisin","Daptomycin",
              "Nonactin","Paraquat","Linezolid","Tigecycline",
              "Spiramycin","Erythromycin","Clarithromycin",
              "Caffeine","Reserpine","Procaine","Theophylline",
              "Verapamil","Meropenem","Triclosan",
              "Gentamicin","Berberine","Metformin","Diclofenac")

mat <- mat |> 
  filter(Drug1 %in% same_res) |> 
  filter(Drug2 %in% same_res)



# TSNE large range
clusters <- read.csv("data/B.SameRes_across_strains_DDIs/3.InteractionScores_tSNE/tSNE_Clustermembership_ppx_large.csv")

mat1 <- full_join(mat,clusters,by="drug_pair") 
mat2 <- mat1[,12:17] %>% t()


mod1 <- phylo.modularity(mat2,partition.gp = mat1$clusters,phy=tree_nw,print.progress=T,
                         CI = TRUE)
mod1
# CR: 1.0139
# P-value: 0.045
# Effect Size: -1.5779
# Based on 1000 random permutations
# Confidence Intervals 0.9857 
# Confidence Intervals 1.0296


# TSNE small range
clusters <- read.csv("data/B.SameRes_across_strains_DDIs/3.InteractionScores_tSNE/tSNE_Clustermembership_ppx_small.csv")

mat1 <- full_join(mat,clusters,by="drug_pair")
mat2 <- mat1[,12:17] %>% t()


mod2 <- phylo.modularity(mat2,partition.gp = mat1$clusters,phy=tree_nw,print.progress=T,
                         CI = TRUE)

mod2

# CR: 0.9897
# P-value: 0.001
# Effect Size: -4.6892
# Based on 1000 random permutations
# Confidence Intervals 0.9195 
# Confidence Intervals 1.0452


# Adams,2017. The most negative effect size (ZCR) is identified, 
# this represented the hypothesis representing the strongest modular signal.

#In this case mod2 is better.

