#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("limma")
#BiocManager::install("DESeq2")

library(DESeq2)
library(scLR)
library(dplyr)
options(scipen = 999)
setwd("C:/Users/farca/Desktop/Project_SingleCellSpatial/Data")

## LR pairs list
lrpairs <- read.table('singlecellsignalR_LRdb.txt', header = TRUE)
lrpairs <- distinct(lrpairs, ligand, receptor, .keep_all= TRUE)

load('COVID19.RData')

#head(counts); dim(counts)
#head(SampleInfo); dim(SampleInfo)
cluster.name <- unique(SampleInfo$cluster)
gene.name <- toupper(rownames(counts))

output <- scLR(countmatrix=counts,
               cellinfo=SampleInfo,
               lrpairs.sample=lrpairs,
               low.filter = 1, p.adjust.method = "BH",
               Brep = 1000, Brep0 = 200, pv0 = 0.05, p.filter = 0.10, show.all = FALSE,
               parallel.use = TRUE, cpucores = 4, zero.impute = TRUE, adjust = 3, rhos = 0,
               normalization = "ByALL", impute.miss.celltype = NA)
RR <- output$Rs


### Only show cases with adj.p < 0.01
FRR <- RR[(RR$adj.p < 0.01),1:11]
#write.csv(FRR, file = "COVID19_ByALL0.csv", row.names = F)


### Show significant lr pairs with means of ligands and receptors > 1
idx.ligand <- cbind(pmatch(sub("(.*)-.*", "\\1",FRR[,2]), gene.name, duplicates.ok = T),
                    pmatch(sub("(.*)-.*", "\\1",FRR[,1]), cluster.name, duplicates.ok = T))
idx.receptor <- cbind(pmatch(sub(".*-", "\\",FRR[,2]), gene.name, duplicates.ok = T),
                    pmatch(sub(".*-", "\\",FRR[,1]), cluster.name, duplicates.ok = T))
aftermean1 <- (output$mu.matrix1[idx.ligand] * 2 + 5 * output$mu.matrix2[idx.ligand])/7
aftermean2 <- (output$mu.matrix1[idx.receptor] * 2 + 5 * output$mu.matrix2[idx.receptor])/7
FRR2 <- FRR[which(aftermean1 > 1 & aftermean2 > 1),1:11]
#write.csv(FRR2, file = "COVID19_ByALL.csv", row.names = F)


### Further exclude Unknown cell type
idx.unknown <- which(sub("(.*)-.*", "\\1",FRR2[,1])=="Unknown" | sub(".*-", "\\",FRR2[,1])=="Unknown")
FRR3 <- FRR2[-idx.unknown,1:11]
write.csv(FRR3, file = "COVID19_ByALL_rmUnknown.csv", row.names = F)




