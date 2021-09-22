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
FRR <- RR[(RR$adj.p < 0.01), ]
#write.csv(FRR, file = "COVID19_ByALL0.csv", row.names = F)


### Show significant lr pairs with means of ligands and receptors > 1
idx.ligand <- cbind(pmatch(sub("(.*)-.*", "\\1",FRR[,2]), gene.name, duplicates.ok = T),
                    pmatch(sub("(.*)-.*", "\\1",FRR[,1]), cluster.name, duplicates.ok = T))
idx.receptor <- cbind(pmatch(sub(".*-", "\\",FRR[,2]), gene.name, duplicates.ok = T),
                    pmatch(sub(".*-", "\\",FRR[,1]), cluster.name, duplicates.ok = T))
aftermean1 <- (output$mu.matrix1[idx.ligand] * 2 + 5 * output$mu.matrix2[idx.ligand])/7
aftermean2 <- (output$mu.matrix1[idx.receptor] * 2 + 5 * output$mu.matrix2[idx.receptor])/7
FRR2 <- FRR[which(aftermean1 > 1 & aftermean2 > 1),]
#write.csv(FRR2, file = "COVID19_ByALL.csv", row.names = F)


### Further exclude Unknown cell type
idx.unknown <- which(sub("(.*)-.*", "\\1",FRR2[,1])=="Unknown" | sub(".*-", "\\",FRR2[,1])=="Unknown")
FRR3 <- FRR2[-idx.unknown,]
write.csv(FRR3, file = "COVID19_ByALL_rmUnknown.csv", row.names = F)



### Figure S10
par(mfrow = c(1,1), mar = c(5,5,3,3))
xx <- -log10(FRR3$adj.limma.Lg.p)
yy <- -log10(FRR3$adj.limma.Rg.p)
plot(xx, yy, pch=16, xlab="-log10FDR (Ligand)", ylab="-log10FDR (Receptor)",
     cex.lab=1.5, cex.axis=1.5)
points(xx[(FRR3$limma.Lg.logFC * FRR3$limma.Rg.logFC) > 0],
       yy[(FRR3$limma.Lg.logFC * FRR3$limma.Rg.logFC) > 0], col="red", pch=16)
abline(h=1, col="blue", lty=2, lwd=2)
abline(v=1, col="blue", lty=2, lwd=2)
FDR <- 0.1
p11 <- round(100*(1 - sum(FRR3$adj.limma.Lg.p<FDR | FRR3$adj.limma.Rg.p<FDR)/1867)) # 16%
p21 <- round(100*(1 - sum(FRR3$adj.limma.Lg.p>FDR | FRR3$adj.limma.Rg.p<FDR)/1867)) # 30%
p12 <- round(100*(1 - sum(FRR3$adj.limma.Lg.p<FDR | FRR3$adj.limma.Rg.p>FDR)/1867)) # 44%
p22 <- round(100*(1 - sum(FRR3$adj.limma.Lg.p>FDR | FRR3$adj.limma.Rg.p>FDR)/1867)) # 10%

text(x=0.2, y=0.2, labels = paste0(p11,"%"), cex=1.5)
text(x=3.3, y=0.2, labels = paste0(p21,"%"), cex=1.5)
text(x=0.2, y=5, labels = paste0(p12,"%"), cex=1.5)
text(x=3, y=5, labels = paste0(p22,"%"), cex=1.5)



