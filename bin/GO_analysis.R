library(data.table)
library(clusterProfiler)
library("org.EcK12.eg.db")

trtab <- fread("results/allddi/idmapping_2023_11_02.tsv",header=TRUE)

clus1 <- fread("results/allddi/targets_cluster.csv",header=TRUE)
clus2 <- lapply(clus1$targets,function(dx)strsplit(dx,split=",")[[1]])
clus4 <- lapply(clus2,function(dx) trtab$RefSeq[match(unlist(dx), trtab$From)])
clus5 <- lapply(clus4,function(dx){
   sapply(dx, function(ex){
    strsplit(ex,split=";")[[1]][1]
      })
    })

clus6 <- lapply(clus5,function(dx){
   sapply(dx, function(ex){
    strsplit(ex,split=".",fixed=TRUE)[[1]][1]
      })
    })
names(clus6) <- c(paste0("cluster",1:16),"all")

x <- enrichGO(
  clus6[[3]],
  OrgDb="org.EcK12.eg.db",
  keyType = "REFSEQ",
  pvalueCutoff = 0.05,
  universe=clus6[[17]],
  qvalueCutoff = 1,
  readable = FALSE,
  pool = FALSE
)


xx <- compareCluster(clus6[1:16],fun="enrichGO",
                     keyType="REFSEQ",
                     OrgDb="org.EcK12.eg.db",
                     universe=clus6[[17]])
pdf("results/allddi/GOfig.pdf",width=8,height = 4)
dotplot(xx, includeAll=TRUE,
        showCategory=50,
        font.size = 8)
dev.off()
