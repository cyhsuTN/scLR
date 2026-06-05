#' Test dysregulated ligand-receptor interactions at small sample sizes (count matrix)
#' @description A function to test dysregulated ligand-receptor interactions from single cell transcriptomics
#' between two conditions.
#' @importFrom Matrix rowSums head tail
#' @importFrom limma lmFit eBayes
#' @importFrom stats rchisq quantile runif t.test sd bw.SJ pnorm p.adjust cor.test wilcox.test
#' @importFrom parallel makeCluster clusterExport parApply stopCluster clusterEvalQ
#' @importFrom dplyr distinct
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors counts
#' @importFrom energy mvnorm.etest

#' @param countmatrix Input. Count matrix/dataframe with gene names.
#' @param cellinfo Input. Information for each cell, including their sampleIDs, conditions, and cell types.
#' @param lrpairs.sample  Input. A K by 2 matrix/dataframe to show names of ligand-receptor pairs
#' which will be compared between two conditions.
#' @param low.filter Input. It is used to pre-exclude LR pairs with lowly expressed.
#' The LR pair will be filtered and not be compared if the unnormalized expression mean of the ligand or
#' the receptor across all samples is less than \emph{low.filter}.
#' The default is 2. If NULL, no LR pairs are pre-excluded.
#' @param p.adjust.method Input (optional). The default is "BH". Other methods include
#  "holm", "hochberg", "hommel", "bonferroni", "BH", "BY" and "fdr" (Please refer to stats::p.adjust).
#' @param Brep  Input (optional). Number of Monte-Carlo simulation replications. The default is 1000.
#' @param Brep0 Input (optional). Number of Monte-Carlo simulation replications in the first-stage (rough) selection.
#' The default is 200. If Brep0 = Brep, no second-stage selection.
#' @param pv0   Input (optional). Prespecified p-value in the first-stage (rough) selection.
#' The default is 0.05.
#' @param p.filter   Input (optional). It is also used to exclude LR pairs with lowly expressed.
#' If the normalized expression means of the ligand or the receptor in both conditions
#' are less than the (p.filter)-quantiles of gene means of the corresponding cell type,
#' the LR pair will be filtered and not be compared. The default is 0.1.
#' @param show.all   Input (optional). A logical indicator whether or not to show
#' the results of the ligand-receptor pairs without being compared. The default is FALSE,
#' showing only the results of the ligand-receptor pairs with being compared.
#' @param parallel.use   Input (optional). A logical indicator whether or not to use a parallel computation.
#' The default is FALSE, i.e., a parallel computation will not be used.
#' @param cpucores   Input (optional). \emph{cpucores} CPU cores (default = 2) are used for parallel computation.
#' @param zero.impute Input (optional). A logical indicator whether zero imputation is performed.
#' If TRUE, when the expressions of ligands or receptors in all samples of one condition are zeros,
#' zeros will be imputed with a constant of 1 before normalization.
#' @param adjust Input (optional). To adjust estimated bandwidth h (bw.SJ * adjust).
#' The default of \emph{adjust} is 3.
#' @param rhos Input (optional). The correlation coefficients in ligand-receptor pairs.
#' The default is 0, assuming zero correlations for all ligand-receptor pairs.
#' Users can assume an overall correlation value (between 0 and 1) for all ligand-receptor pairs, e.g, rhos = 0.1;
#' or set rhos = "est", the Pearson correlation coefficient estimates are used for the significant correlations
#' (p-values < 0.01 determined by R function cor.test), and 0 for the non-significant correlations.
#' @param normalization Input (optional). Perform normalization within each cell type (normalization = "ByCT") or
#' for all cell types together (normalization = "ByALL"). The default is "ByALL".
#' @param impute.miss.celltype Input (optional). Sometimes, some samples do not have specific cell types.
#' For the missing cell type in the samples, the gene expressions will be replaced with \emph{impute.miss.celltype}.
#' The default is NA, denoting there not exist the cell types; 0, denoting lowly expressed.
#' @param bntest Input (optional). A logical indicator whether or not to perform a bivariate noraml test.
#' The default is FALSE.
#' @param sig.bntest Input (optional). If bntest = TRUE, under what significance level the significant pairs will be further tested
#' by a bivariate noraml test. The default is 0.05.
#'
#'

#' @return \item{obs.xy.diff}{Observed difference of mean of LR between two conditions (condition 2 minus condition 1)}
#' @return \item{null.diff.sd}{Standard deviation of mean difference of LR under null hypothesis}
#' @return \item{diff.stat}{obs.xy.diff is divided by null.diff.sd}
#' @return \item{pvalue}{P-value}
#' @return \item{stage}{Which stage gene pairs pass to}
#' @return \item{adj.p}{Adjusted p-value of scLR for multiple comparison}
#' @return \item{Welch.t.*}{Statistic, pooled standard error, and p-value of Welch's t-test}
#' @return \item{Welch.t.adj.p}{Adjusted p-value of Welch's t-test for multiple comparison}
#' @return \item{WRS.*}{Statistic and p-value of Wilcoxon Rank Sum test}
#' @return \item{WRS.adj.p}{Adjusted p-value of Wilcoxon Rank Sum test for multiple comparison}
#' @return \item{bntest.p}{p-value of bivariate normal test if bntest = TRUE}
#' @return \item{bntest.adj.p}{Adjusted p-value of bivariate normal test for multiple comparison if bntest = TRUE}
#' @return \item{mu.matrix1}{Sample mean for each gene in condition 1}
#' @return \item{mu.matrix2}{Sample mean for each gene in condition 2}
#' @return \item{sigma.matrix1}{Posterior residual standard deviation for each gene in condition 1}
#' @return \item{sigma.matrix2}{Posterior residual standard deviation for each gene in condition 2}


#' @examples # A simulated data to compare LR pairs between two conditions: TX1 (control) and TX2.
#' @examples # 3 replicates each condition.
#' @examples # 1000 genes and 100 cells (5 cell types, 20 cells each) each replicate.
#' @examples # 10 ligand-receptor gene pairs across 5 cell types are compared between TX1 and TX2.
#'
#' @examples set.seed(2022)
#' @examples cellinfo <- data.frame(sampleID = factor(paste0("s", rep(1:6, each=100))),
#' @examples                       condition = factor(paste0("tx", rep(1:2, each=300))),
#' @examples                       cellcluster = factor(paste0("cc", rep(rep(1:5, each=20), 6))) )
#' @examples cellinfo$condition <- relevel(cellinfo$condition, ref="tx1") # Let tx1 be control.
#' @examples G <- 1000; n <- 600 # 1000 genes, 600 cells
#' @examples # Data are generated from NB dist.
#' @examples mu1 <- rgamma(G, shape = 2, rate = 2)
#' @examples NB_cell <- function(j) rnbinom(G, size = 0.1, mu = mu1)
#' @examples countmatrix <- as(sapply(1:n, NB_cell),"sparseMatrix"); rownames(countmatrix) <- paste0("g", 1:1000)
#' @examples # Names of 10 ligand-receptor pairs
#' @examples lrpairs.sample <- data.frame(ligand=paste0("g", sample(1:1000, 10)),
#' @examples                             receptor=paste0("g", sample(1:1000, 10)) )
#' @examples output <- scLR(countmatrix, cellinfo, lrpairs.sample)
#' @examples head(output$Rs, 20)
#'
#' @examples # Assume sample s2 does NOT have cell type cc3
#' @examples idx.remove <- which(cellinfo$sampleID=="s2" & cellinfo$cellcluster=="cc3")
#' @examples cellinfo1 <- cellinfo[-idx.remove,]; cellinfo1 <- droplevels(cellinfo1)
#' @examples countmatrix1 <- countmatrix[,-idx.remove]
#' @examples output <- scLR(countmatrix1, cellinfo1, lrpairs.sample, impute.miss.celltype = NA)
#' @examples head(output$Rs, 20)

#' @export

scLR <- function(countmatrix, cellinfo, lrpairs.sample, low.filter = 1, p.adjust.method = "BH",
                  Brep = 1000, Brep0 = 200, pv0 = 0.05, p.filter = 0.1, show.all = FALSE,
                  parallel.use = FALSE, cpucores = 2, zero.impute = TRUE, adjust = 3, rhos = 0,
                  normalization = "ByALL", impute.miss.celltype = NA, bntest = FALSE, sig.bntest = 0.01) {


  colnames(cellinfo) <- c("sampleID", "condition", "cellcluster")
  colnames(lrpairs.sample) <- c("ligand", "receptor")
  rownames(countmatrix) <- toupper(rownames(countmatrix))
  lrpairs.sample$ligand <- toupper(lrpairs.sample$ligand)
  lrpairs.sample$receptor <- toupper(lrpairs.sample$receptor)

  ### logic check
  if(ncol(countmatrix) != nrow(cellinfo)) {
    stop("ncol(countmatrix) != nrow(cellinfo)")
  }

  ### sum counts of same cell cluster by samples
  samples <- unique(cellinfo$sampleID)
  cellclusters <- unique(cellinfo$cellcluster)

  nsamples <- length(samples)
  No_celltypes <- length(cellclusters)
  no.genes <- nrow(countmatrix)
  if(nsamples < 2) {
    stop("Only one sample in this dataset")
  }

  conditions <- unlist(lapply(samples, function(ss) unique(cellinfo$condition[cellinfo$sampleID==ss])))

  sum.same.cc <- function(sampleid, ccname) {
    idx <- which(cellinfo$sampleID==sampleid & cellinfo$cellcluster==ccname)
    if (length(idx)==0) {
      rep(impute.miss.celltype, no.genes)
    } else {
      rowSums(countmatrix[,which(cellinfo$sampleID==sampleid & cellinfo$cellcluster==ccname), drop=FALSE])
    }
  }

  y <- list(sumExp=NULL, group=NULL)
  for(sampleid in samples) {
    aaa <- sapply(cellclusters, function(cc) sum.same.cc(sampleid, cc) )
    y$sumExp <- cbind(y$sumExp, aaa)
  }
  y$group <- data.frame(sampleID=rep(samples, each=No_celltypes),
                        condition=rep(conditions, each=No_celltypes),
                        cellcluster=rep(cellclusters, nsamples))

  ### remove data which will not be used in subsequent computation
  rm(countmatrix, cellinfo); gc()



  ### normalize counts and take log2( + 1)
  signal.exp <- y$sumExp
  if (normalization != "ByALL") {
    for(i in 1:No_celltypes) {
      cell.idx0 <- (0:(nsamples-1))*No_celltypes + i
      value.colsum <- colSums(y$sumExp[,cell.idx0])
      idx.na <- which(is.na(value.colsum) | value.colsum==0)
      if (length(idx.na)>0) {
        cell.idx <- cell.idx0[-idx.na]
        dds.c <- DESeqDataSetFromMatrix(countData = y$sumExp[,cell.idx],
                                        colData = y$group[cell.idx,],
                                        design = ~ condition)
        dds.c <- estimateSizeFactors(dds.c)
        signal.exp[,cell.idx] <- log2(counts(dds.c, normalized = TRUE)+1)
      } else {
        cell.idx <- cell.idx0
        dds.c <- DESeqDataSetFromMatrix(countData = y$sumExp[,cell.idx],
                                        colData = y$group[cell.idx,],
                                        design = ~ condition)
        dds.c <- estimateSizeFactors(dds.c)
        signal.exp[,cell.idx] <- log2(counts(dds.c, normalized = TRUE)+1)
      }
    }
  } else {
    cell.idx0 <- 1:ncol(signal.exp)
    value.colsum <- colSums(y$sumExp[,cell.idx0])
    idx.na <- which(is.na(value.colsum) | value.colsum==0)
    if (length(idx.na)>0) {
      cell.idx <- cell.idx0[-idx.na]
      dds.c <- DESeqDataSetFromMatrix(countData = y$sumExp[,cell.idx],
                                      colData = y$group[cell.idx,],
                                      design = ~ condition)
      dds.c <- estimateSizeFactors(dds.c)
      signal.exp[,cell.idx] <- log2(counts(dds.c, normalized = TRUE)+1)
    } else {
      cell.idx <- cell.idx0
      dds.c <- DESeqDataSetFromMatrix(countData = y$sumExp[,cell.idx],
                                      colData = y$group[cell.idx,],
                                      design = ~ condition)
      dds.c <- estimateSizeFactors(dds.c)
      signal.exp[,cell.idx] <- log2(counts(dds.c, normalized = TRUE)+1)
    }
  }

  ### Two sample matrices
  if(is.null(levels(y$group$condition))) {
    condition.type <- unique(y$group$condition)
  } else {
    condition.type <- levels(y$group$condition)
  }
  samplematrix1 <- signal.exp[,which(y$group$condition==condition.type[1])]
  samplematrix2 <- signal.exp[,which(y$group$condition==condition.type[2])]
  colnames(samplematrix1) <- y$group$cellcluster[which(y$group$condition==condition.type[1])]
  colnames(samplematrix2) <- y$group$cellcluster[which(y$group$condition==condition.type[2])]

  ### remove lr pairs which can not be found in two sample matrices
  idx.genename.nofound <- which(!(lrpairs.sample$ligand %in% rownames(samplematrix1)) |
                                  !(lrpairs.sample$receptor %in% rownames(samplematrix1)))
  if(sum(idx.genename.nofound)==0) {
    lrpairs.v2 <- lrpairs.sample
  } else {
    lrpairs.v2 <- lrpairs.sample[-c(idx.genename.nofound),]
  }
  genematrix <- data.frame(pmatch(lrpairs.v2$ligand, rownames(samplematrix1), duplicates.ok = T),
                           pmatch(lrpairs.v2$receptor, rownames(samplematrix1), duplicates.ok = T))
  genematrix <- distinct(genematrix, .keep_all= TRUE)


  ### define a genepair matrix
  genepairs <- NULL
  for (i in 1:No_celltypes) {
    for (j in 1:No_celltypes) {
      genepairs <- rbind(genepairs, cbind(i, j, genematrix))
    }
  }
  genepairs <- genepairs[,c(1,3,2,4)]


  ### check if the unnormalized expression means of ligands or receptors across all samples
  ### are less than low.filter (before normalization)
  ### If yes, the LR pairs will not be compared
  if(!is.null(low.filter) & low.filter > min(y$sumExp, na.rm=TRUE)) {
    genes.exclude <- NULL
    for(i in 1:No_celltypes) {
      cell.idx <- (0:(nsamples-1))*No_celltypes + i
      idx.small <- which(rowMeans(y$sumExp[,cell.idx], na.rm=T) < low.filter)
      genes.exclude <- rbind(genes.exclude, cbind(rep(i,length(idx.small)), idx.small))
    }
    if(nrow(genes.exclude)>0) {
      genes.exclude.connect <- paste0(genes.exclude[,1],"-",genes.exclude[,2])
      genepairs.l.connect <- paste0(genepairs[,1],"-",genepairs[,2])
      genepairs.r.connect <- paste0(genepairs[,3],"-",genepairs[,4])

      idx.genename.exclude <- which( (genepairs.l.connect %in% genes.exclude.connect) |
                                       (genepairs.r.connect %in% genes.exclude.connect)  )
      if(length(idx.genename.exclude)>0) {
        genepairs <- genepairs[-c(idx.genename.exclude),]
      }
    }
  }

  #set.seed(2021)
  RR.all <- scXY(samplematrix1, samplematrix2, genepairs, No_celltypes,
                 Brep, Brep0, pv0, p.filter, show.all,
                 parallel.use, cpucores, zero.impute, adjust, rhos)

  RR.all$Rs$adj.p <- round(p.adjust(as.numeric(as.character(RR.all$Rs$pvalue)), method=p.adjust.method), 12)

  RR.all$Rs$Welch.t.p[which(is.na(RR.all$Rs$Welch.t.p))] <- 1
  RR.all$Rs$Welch.t.adj.p <- round(p.adjust(as.numeric(as.character(RR.all$Rs$Welch.t.p)), method=p.adjust.method), 12)

  RR.all$Rs$WRS.p[which(is.na(RR.all$Rs$WRS.p))] <- 1
  RR.all$Rs$WRS.adj.p <- round(p.adjust(as.numeric(as.character(RR.all$Rs$WRS.p)), method=p.adjust.method), 12)


  #### calculate adjusted p-values of using limma approach
  idx.for.limma <- as.matrix(
    distinct(
      data.frame(rbind(as.matrix(RR.all$Rs[,c(2,1)]), as.matrix(RR.all$Rs[, c(4,3)]))),
      .keep_all = TRUE)
  )
  colnames(idx.for.limma) <- c("gene", "ctype")

  adj.limma.p <- round(matrix(p.adjust(as.numeric(RR.all$pv.matrix), method=p.adjust.method),
                              nrow=nrow(RR.all$pv.matrix))[idx.for.limma], 12)
    #adj.limma.p <- round(p.adjust(as.numeric(RR.all$pv.matrix[idx.for.limma]), method=p.adjust.method), 12)
  gc.connect <- paste0(idx.for.limma[,"gene"],"-",idx.for.limma[,"ctype"])

  RR.all$Rs$limma.Lg.adj.p <- adj.limma.p[pmatch(paste0(RR.all$Rs$ligand_gene,"-",RR.all$Rs$ligand_cell), gc.connect, duplicates.ok = T)]
  RR.all$Rs$limma.Rg.adj.p <- adj.limma.p[pmatch(paste0(RR.all$Rs$receptor_gene,"-",RR.all$Rs$receptor_cell), gc.connect, duplicates.ok = T)]


  sample.group <- lapply(condition.type, function(x) as.character(unique(y$group$sampleID[which(y$group$condition==x)])))
  names(sample.group) <- condition.type

  #sig.bntest <- 0.95
  ## bivariate normal test  # library(energy)
  if (bntest) {
    idx.bntest <- which(RR.all$Rs$adj.p < sig.bntest)
    if(length(idx.bntest) > 0) {

      RR.all$Rs$bntest.p <- RR.all$Rs$bntest.adj.p <- NA

      N1 <- ncol(samplematrix1)/No_celltypes
      N2 <- ncol(samplematrix2)/No_celltypes
      const1 <- (0:(N1 - 1)) * No_celltypes
      const2 <- (0:(N2 - 1)) * No_celltypes

      sig.pairs.bntest <- RR.all$Rs[idx.bntest, 1:4]
      bntest.p <- apply(sig.pairs.bntest, 1, function(pair) {
        dat1 <- data.frame(L = c(samplematrix1[pair[2], pair[1]+const1]), R = c(samplematrix1[pair[4], pair[3]+const1]))
        dat2 <- data.frame(L = c(samplematrix2[pair[2], pair[1]+const2]), R = c(samplematrix2[pair[4], pair[3]+const2]))

        bntest1 <- mvnorm.etest(dat1[!is.na(dat1$L) & !is.na(dat1$R),], R=100)
        bntest2 <- mvnorm.etest(dat2[!is.na(dat2$L) & !is.na(dat2$R),], R=100)
        min(2 * min(bntest1$p.value, bntest2$p.value, na.rm = T), 1, na.rm = T)
      })
      RR.all$Rs$bntest.p[idx.bntest] <- bntest.p
      RR.all$Rs$bntest.adj.p[idx.bntest] <- round(p.adjust(bntest.p, method=p.adjust.method), 12)

      bntest.data <- data.frame(cbind(bntest.p, bntest.adj.p = round(p.adjust(bntest.p, method=p.adjust.method), 12),
                                      t(apply(sig.pairs.bntest, 1, function(pair) {
                                        c(L = c(samplematrix1[pair[2], pair[1]+const1], samplematrix2[pair[2], pair[1]+const2]),
                                          R = c(samplematrix1[pair[4], pair[3]+const1], samplematrix2[pair[4], pair[3]+const2]))
                                      }))))
      colnames(bntest.data)[-c(1:2)] <- c(paste0(unlist(sample.group), c(".L")), paste0(unlist(sample.group), c(".R")))

      lr.cell.name <- paste0(cellclusters[sig.pairs.bntest[,1]],"-",cellclusters[sig.pairs.bntest[,3]])
      lr.gene.name <- paste0(rownames(samplematrix1)[sig.pairs.bntest[,2]],"-",rownames(samplematrix1)[sig.pairs.bntest[,4]])
      bntest.data <- cbind(lr.cell.name, lr.gene.name, bntest.data)
    } else {
      bntest.data <- NULL
    }
  } else {
      bntest.data <- NULL
  }

  N1 <- ncol(samplematrix1)/No_celltypes
  N2 <- ncol(samplematrix2)/No_celltypes
  const1 <- (0:(N1 - 1)) * No_celltypes
  const2 <- (0:(N2 - 1)) * No_celltypes
  sig.pairs.bntest <- RR.all$Rs[, 1:4]
  sumRawData <- data.frame(cbind(#NA, bntest.adj.p = NA,
                                  t(apply(sig.pairs.bntest, 1, function(pair) {
                                    c(L = c(samplematrix1[pair[2], pair[1]+const1], samplematrix2[pair[2], pair[1]+const2]),
                                      R = c(samplematrix1[pair[4], pair[3]+const1], samplematrix2[pair[4], pair[3]+const2]))
                                  }))))
  colnames(sumRawData) <- c(paste0(unlist(sample.group), c(".L")), paste0(unlist(sample.group), c(".R")))



  lr.cell.name <- paste0(cellclusters[RR.all$Rs$ligand_cell],"-",cellclusters[RR.all$Rs$receptor_cell])
  lr.gene.name <- paste0(rownames(samplematrix1)[RR.all$Rs$ligand_gene],"-",rownames(samplematrix1)[RR.all$Rs$receptor_gene])


  if(bntest & !is.null(bntest.data)) {
    FRR2 <- cbind(lr.cell.name, lr.gene.name,
                  RR.all$Rs[,c("obs.xy.diff", "null.diff.sd", "pvalue", "stage", "adj.p",
                               "Welch.t.stat", "Welch.t.sd", "Welch.t.p", "Welch.t.adj.p",
                               "WRS.p", "WRS.adj.p", "bntest.p", "bntest.adj.p",
                               "limma.Lg.p", "limma.Rg.p", "limma.Lg.logFC", "limma.Rg.logFC",
                               "limma.Lg.adj.p", "limma.Rg.adj.p")
                            ])
  } else {
    FRR2 <- cbind(lr.cell.name, lr.gene.name, RR.all$Rs[,c("obs.xy.diff", "null.diff.sd", "pvalue", "stage", "adj.p",
                                                           "Welch.t.stat", "Welch.t.sd", "Welch.t.p", "Welch.t.adj.p",
                                                           "WRS.p", "WRS.adj.p",
                                                           "limma.Lg.p", "limma.Rg.p", "limma.Lg.logFC", "limma.Rg.logFC",
                                                           "limma.Lg.adj.p", "limma.Rg.adj.p")])
  }


  list(Rs=FRR2,
       conditions = paste0("Condition_",1:2,": ",condition.type),
       sample.group = sample.group,
       bntest.data = bntest.data,
       sumRawData = sumRawData,
       mu.matrix1=RR.all$mu.matrix1,
       mu.matrix2=RR.all$mu.matrix2,
       sigma.matrix1=RR.all$sigma.matrix1,
       sigma.matrix2=RR.all$sigma.matrix2,
       pv.matrix=RR.all$pv.matrix,
       logFC.matrix=RR.all$logFC.matrix
  )

}

#RR.all <- output1
#relative.difference(RR.all)
relative.difference <- function(RR.all) {

  a1 <- ((RR.all$mu.matrix1[as.matrix(RR.all$Rs[,c("ligand_gene", "ligand_cell")])]/
            RR.all$sigma.matrixNull[as.matrix(RR.all$Rs[,c("ligand_gene", "ligand_cell")])] +
            RR.all$mu.matrix1[as.matrix(RR.all$Rs[,c("receptor_gene", "receptor_cell")])]/
            RR.all$sigma.matrixNull[as.matrix(RR.all$Rs[,c("receptor_gene", "receptor_cell")])]) / sqrt(2*(1+RR.all$Rs$rhoNull)) )^2

  b1 <- ((RR.all$mu.matrix1[as.matrix(RR.all$Rs[,c("ligand_gene", "ligand_cell")])]/
            RR.all$sigma.matrixNull[as.matrix(RR.all$Rs[,c("ligand_gene", "ligand_cell")])] -
            RR.all$mu.matrix1[as.matrix(RR.all$Rs[,c("receptor_gene", "receptor_cell")])]/
            RR.all$sigma.matrixNull[as.matrix(RR.all$Rs[,c("receptor_gene", "receptor_cell")])]) / sqrt(2*(1-RR.all$Rs$rhoNull)) )^2

  a2 <- ((RR.all$mu.matrix2[as.matrix(RR.all$Rs[,c("ligand_gene", "ligand_cell")])]/
            RR.all$sigma.matrixNull[as.matrix(RR.all$Rs[,c("ligand_gene", "ligand_cell")])] +
            RR.all$mu.matrix2[as.matrix(RR.all$Rs[,c("receptor_gene", "receptor_cell")])]/
            RR.all$sigma.matrixNull[as.matrix(RR.all$Rs[,c("receptor_gene", "receptor_cell")])]) / sqrt(2*(1+RR.all$Rs$rhoNull)) )^2

  b2 <- ((RR.all$mu.matrix2[as.matrix(RR.all$Rs[,c("ligand_gene", "ligand_cell")])]/
            RR.all$sigma.matrixNull[as.matrix(RR.all$Rs[,c("ligand_gene", "ligand_cell")])] -
            RR.all$mu.matrix2[as.matrix(RR.all$Rs[,c("receptor_gene", "receptor_cell")])]/
            RR.all$sigma.matrixNull[as.matrix(RR.all$Rs[,c("receptor_gene", "receptor_cell")])]) / sqrt(2*(1-RR.all$Rs$rhoNull)) )^2

  sigmaLR1 <- RR.all$sigma.matrix1[as.matrix(RR.all$Rs[,c("ligand_gene", "ligand_cell")])] *
    RR.all$sigma.matrix1[as.matrix(RR.all$Rs[,c("receptor_gene", "receptor_cell")])]

  sigmaLR2 <- RR.all$sigma.matrix2[as.matrix(RR.all$Rs[,c("ligand_gene", "ligand_cell")])] *
    RR.all$sigma.matrix2[as.matrix(RR.all$Rs[,c("receptor_gene", "receptor_cell")])]

  muLR1 <- RR.all$mu.matrix1[as.matrix(RR.all$Rs[,c("ligand_gene", "ligand_cell")])] *
    RR.all$mu.matrix1[as.matrix(RR.all$Rs[,c("receptor_gene", "receptor_cell")])]

  muLR2 <- RR.all$mu.matrix2[as.matrix(RR.all$Rs[,c("ligand_gene", "ligand_cell")])] *
    RR.all$mu.matrix2[as.matrix(RR.all$Rs[,c("receptor_gene", "receptor_cell")])]

  data.frame(#RD.muLR = (muLR2 - muLR1)/((muLR1 + muLR2)/2),
        RD.a = (a2 - a1)/((a1 + a2)/2),
        RD.b = (b2 - b1)/((b1 + b2)/2),
        RD.sigmaLR = (sigmaLR2 - sigmaLR1)/((sigmaLR1 + sigmaLR2)/2))


}



scXY <- function(samplematrix1, samplematrix2, genepairs, No_celltypes,
                 Brep = 1000, Brep0 = 200, pv0 = 0.05, p.filter = 0.1, show.all = FALSE,
                 parallel.use = FALSE, cpucores = 2, zero.impute = TRUE, adjust = 3, rhos = 0) {

  if(!is.matrix(genepairs)) {genepairs <- as.matrix(genepairs)}

  N1 <- ncol(samplematrix1)/No_celltypes
  N2 <- ncol(samplematrix2)/No_celltypes

  ### pairname convert to index
  genepairs2 <- genepairs

  mu.matrix1 <- sigma.matrixNull <- sigma.matrix1D <- matrix(NA, nrow = nrow(samplematrix1), ncol = No_celltypes)
  mu.matrix2 <- sigma.matrix2D <- matrix(NA, nrow = nrow(samplematrix2), ncol = No_celltypes)
  pv.matrix <- logFC.matrix <- matrix(NA, nrow = nrow(samplematrix1), ncol = No_celltypes)
  const1 <- (0:(N1 - 1)) * No_celltypes
  const2 <- (0:(N2 - 1)) * No_celltypes

  min1 <- apply(samplematrix1, 2, function(x) { ifelse(all(is.na(x)), NA, min(x[x>0], na.rm = T)/4) })
  min2 <- apply(samplematrix2, 2, function(x) { ifelse(all(is.na(x)), NA, min(x[x>0], na.rm = T)/4) })

  ### run limma to estimate mean and posterior residual standard deviations
  set.seed(2021)
  for (k in 1:No_celltypes) {

    idx.col1 <- k + const1
    idx.col2 <- k + const2

    if(zero.impute) {
      idx.allzero1 <- which(rowSums(samplematrix1[, idx.col1, drop=FALSE], na.rm = T)==0)
      idx.allzero2 <- which(rowSums(samplematrix2[, idx.col2, drop=FALSE], na.rm = T)==0)
      min1s <- rep(min1[idx.col1], each=length(idx.allzero1))
      min2s <- rep(min2[idx.col2], each=length(idx.allzero2))
      samplematrix1[idx.allzero1, idx.col1] <- ifelse(is.na(min1s), NA, suppressWarnings(runif(length(min1s), 0, max=min1s)))
      samplematrix2[idx.allzero2, idx.col2] <- ifelse(is.na(min2s), NA, suppressWarnings(runif(length(min2s), 0, max=min2s)))
    }

    ###### allow one sample in one of two conditions
    yy <- cbind(samplematrix1[, idx.col1], samplematrix2[, idx.col2])
    designmatrix <- cbind(Grp1=1, Grp2vs1=c(rep(0, N1), rep(1, N2)))

    fit0 <- lmFit(yy, designmatrix)

    if(fit0$df.residual[1]!=0) {
      fit0 <- eBayes(fit0)
      sigma.matrixNull[, k] <- sqrt(fit0$s2.post)
      pv.matrix[, k]    <-  as.numeric(fit0$p.value[,2])
    } else {
      sigma.matrixNull[, k] <- NA
      pv.matrix[, k]    <-  NA
    }

    mu.matrix1[, k] <- as.numeric(fit0$coefficients[,1])
    mu.matrix2[, k] <- as.numeric(rowSums(fit0$coefficients))
    logFC.matrix[, k] <-  as.numeric(fit0$coefficients[,2])

    ######
    fit1 <- lmFit(samplematrix1[, idx.col1])
    fit2 <- lmFit(samplematrix2[, idx.col2])

    fit1$sigma[which(fit1$sigma==0)] <- 1E-8
    fit2$sigma[which(fit2$sigma==0)] <- 1E-8

    if (fit1$df.residual[1]!=0) {
      fit1 <- eBayes(fit1)
      sigma.matrix1D[, k] <- sqrt(fit1$s2.post)
    } else {
      sigma.matrix1D[, k] <- NA
    }

    if (fit2$df.residual[1]!=0) {
      fit2 <- eBayes(fit2)
      sigma.matrix2D[, k] <- sqrt(fit2$s2.post)
    } else {
      sigma.matrix2D[, k] <- NA
    }

  }

  ### Obtain the thresholds for each cell type in two conditions (after normalization)
  threshold1 <- apply(mu.matrix1, 2, quantile, probs = p.filter)
  threshold2 <- apply(mu.matrix2, 2, quantile, probs = p.filter)

  comparefun <- function(pair) {
    #pair <- genepairs[284,]
    XY1 <- samplematrix1[pair[2], pair[1]+const1] * samplematrix1[pair[4], pair[3]+const1]
    XY2 <- samplematrix2[pair[2], pair[1]+const2] * samplematrix2[pair[4], pair[3]+const2]
    n1 <- sum(!is.na(XY1))
    n2 <- sum(!is.na(XY2))

    ## Check if the normalized expression means of ligands or receptors in both conditions are less than pre-specified quantiles.
    cond1 <- (mu.matrix1[pair[2], pair[1]] <= threshold1[pair[1]]) & (mu.matrix2[pair[2], pair[1]] <= threshold2[pair[1]])
    cond2 <- (mu.matrix1[pair[4], pair[3]] <= threshold1[pair[3]]) & (mu.matrix2[pair[4], pair[3]] <= threshold2[pair[3]])

    if (n1 == 0 | n2 == 0) {
      c(obs.xy.diff = "Only one condition has non-NA XY", null.diff.sd = "-", diff.stat = "-",
        pvalue = "-", stage = "-", rhoNull = "-",
        Welch.t.stat = "-", Welch.t.sd = "-", Welch.t.p = "-", WRS.stat = "-", WRS.p = "-",
        limma.Lg.p = "-", limma.Rg.p = "-", limma.Lg.logFC = "-", limma.Rg.logFC = "-")
    } else if (any(is.na(c(sigma.matrixNull[pair[2], pair[1]], sigma.matrixNull[pair[4], pair[3]])))) {
      c(obs.xy.diff = "NA sigma", null.diff.sd = "-", diff.stat = "-",
        pvalue = "-", stage = "-", rhoNull = "-",
        Welch.t.stat = "-", Welch.t.sd = "-", Welch.t.p = "-", WRS.stat = "-", WRS.p = "-",
        limma.Lg.p = "-", limma.Rg.p = "-", limma.Lg.logFC = "-", limma.Rg.logFC = "-")
    } else if (cond1 | cond2) {
      c(obs.xy.diff = paste0("Lowly expressed"), null.diff.sd = "-", diff.stat = "-",
        pvalue = "-", stage = "-", rhoNull = "-",
        Welch.t.stat = "-", Welch.t.sd = "-", Welch.t.p = "-", WRS.stat = "-", WRS.p = "-",
        limma.Lg.p = "-", limma.Rg.p = "-", limma.Lg.logFC = "-", limma.Rg.logFC = "-")
    } else {
      w1 <- n1/(n1 + n2)
      w2 <- 1 - w1
      obs.xy.diff <- mean(XY2, na.rm = T) - mean(XY1, na.rm = T)

      if (n1 > 1 & n2 > 1 & (sd(XY1, na.rm = T) > 0 | sd(XY2, na.rm = T) > 0) ) {
        Welch.t <- t.test(XY2, XY1) ## Welch t test
        Welch.t$Welch.t.stat <- round(as.numeric(Welch.t$statistic),4)
        Welch.t$Welch.t.sd <- round(as.numeric(Welch.t$stderr),4)
        Welch.t$Welch.t.p <- round(as.numeric(Welch.t$p.value),6)

        WRS.test <- wilcox.test(XY2, XY1) ## Wilcoxon Rank Sum
        WRS.test$WRS.stat <- round(as.numeric(WRS.test$statistic),4)
        WRS.test$WRS.p <- round(as.numeric(WRS.test$p.value),6)
      } else {
        Welch.t <- NULL
        Welch.t$Welch.t.stat <- Welch.t$Welch.t.sd <- Welch.t$Welch.t.p <- NA

        WRS.test <- NULL
        WRS.test$WRS.stat <- WRS.test$WRS.p <- NA
      }

      mu1 <- c(mu.matrix1[pair[2], pair[1]], mu.matrix1[pair[4], pair[3]])
      mu2 <- c(mu.matrix2[pair[2], pair[1]], mu.matrix2[pair[4], pair[3]])

      sigmaNull <- c(sigma.matrixNull[pair[2], pair[1]], sigma.matrixNull[pair[4], pair[3]])
      #sigma1 <- c(sigma.matrix1[pair[2], pair[1]], sigma.matrix1[pair[4], pair[3]])
      #sigma2 <- c(sigma.matrix2[pair[2], pair[1]], sigma.matrix2[pair[4], pair[3]])

      sigma.prodNull <- prod(sigmaNull)
      #sigma.prod1 <- prod(sigma1)
      #sigma.prod2 <- prod(sigma2)

      if(rhos==0) {
        rho1 <- rho2 <- 0
      } else if (rhos=="est") {
        rhoest <- cor.test(c(samplematrix1[pair[2], pair[1]+const1], samplematrix2[pair[2], pair[1]+const2]),
                           c(samplematrix1[pair[4], pair[3]+const1], samplematrix2[pair[4], pair[3]+const2]),
                           method = "pearson")
        if(rhoest$p.value < 0.01) {
          rho1 <- rho2 <- as.numeric(rhoest$estimate)
        } else {
          rho1 <- rho2 <- 0
        }
      } else {
        rho1 <- rho2 <- rhos
      }

      #a1D <- sum(mu1/sigma1)^2/(2*(1+rho1))
      #b1D <- diff(mu1/sigma1)^2/(2*(1-rho1))
      #a2D <- sum(mu2/sigma2)^2/(2*(1+rho2))
      #b2D <- diff(mu2/sigma2)^2/(2*(1-rho2))

      a1 <- sum(mu1/sigmaNull)^2/(2*(1+rho1))
      b1 <- diff(mu1/sigmaNull)^2/(2*(1-rho1))
      a2 <- sum(mu2/sigmaNull)^2/(2*(1+rho2))
      b2 <- diff(mu2/sigmaNull)^2/(2*(1-rho2))
      a.bar <- (w1*a1 + w2*a2)
      b.bar <- (w1*b1 + w2*b2)
      sigma.prod.bar <- sigma.prodNull
      rho.bar <- (w1*rho1 + w2*rho2)

      if (rho1==1) {
        xy.mean.diff.null1 <- function(n1, n2, sigma.prod.bar, a.bar, b.bar, rho.bar) {
          xy <- sigma.prod.bar*rchisq(n1+n2, df=1, ncp=a.bar)
          mean(tail(xy,n2)) - mean(head(xy,n1))
        }
      } else if (rho1==-1) {
        xy.mean.diff.null1 <- function(n1, n2, sigma.prod.bar, a.bar, b.bar, rho.bar) {
          xy <- sigma.prod.bar*(-rchisq(n1+n2, df=1, ncp=b.bar))
          mean(tail(xy,n2)) - mean(head(xy,n1))
        }
      } else {
        xy.mean.diff.null1 <- function(n1, n2, sigma.prod.bar, a.bar, b.bar, rho.bar) {
          xy <- 0.25*sigma.prod.bar*(2*(1+rho.bar)*rchisq(n1+n2, df=1, ncp=a.bar) -
                                       2*(1-rho.bar)*rchisq(n1+n2, df=1, ncp=b.bar))
          mean(tail(xy,n2)) - mean(head(xy,n1))
        }
      }

      #set.seed(2021)
      Ts0 <- replicate(Brep0, xy.mean.diff.null1(n1, n2, sigma.prod.bar, a.bar, b.bar, rho.bar))
      pvalue0 <- mean(abs(Ts0) > abs(obs.xy.diff), na.rm=T)
      if (pvalue0 < pv0 & Brep0 < Brep) {
        Ts <- replicate(Brep-Brep0, xy.mean.diff.null1(n1, n2, sigma.prod.bar, a.bar, b.bar, rho.bar))

        bw <- bw.SJ(c(Ts0,Ts), nb=1000, method = c("ste")) * adjust
        tailprob <- mean(pnorm( obs.xy.diff - c(Ts0,Ts), 0, bw))
        pvalue <- ifelse(tailprob < 0.5, 2 * tailprob, 2 * (1 - tailprob))
        null.diff.sd <- sd(c(Ts0,Ts))
        diff.stat <- obs.xy.diff/null.diff.sd
        stage <- 2
        c(obs.xy.diff = round(obs.xy.diff,4), null.diff.sd = round(null.diff.sd,4), diff.stat = round(diff.stat,4),
          pvalue = round(pvalue,12), stage = stage, rhoNull = rho.bar,
          Welch.t.stat = Welch.t$Welch.t.stat, Welch.t.sd = Welch.t$Welch.t.sd, Welch.t.p = Welch.t$Welch.t.p,
          WRS.stat = WRS.test$WRS.stat, WRS.p = WRS.test$WRS.p,
          limma.Lg.p = round(pv.matrix[pair[2],pair[1]], 12),
          limma.Rg.p = round(pv.matrix[pair[4],pair[3]], 12),
          limma.Lg.logFC = round(logFC.matrix[pair[2],pair[1]],4),
          limma.Rg.logFC = round(logFC.matrix[pair[4],pair[3]],4))
      } else {
        null.diff.sd <- sd(Ts0)
        diff.stat <- obs.xy.diff/null.diff.sd
        stage <- 1
        c(obs.xy.diff = round(obs.xy.diff,4), null.diff.sd = round(null.diff.sd,4), diff.stat = round(diff.stat,4),
          pvalue = pvalue0, stage = stage, rhoNull = rho.bar,
          Welch.t.stat = Welch.t$Welch.t.stat, Welch.t.sd = Welch.t$Welch.t.sd, Welch.t.p = Welch.t$Welch.t.p,
          WRS.stat = WRS.test$WRS.stat, WRS.p = WRS.test$WRS.p,
          limma.Lg.p = round(pv.matrix[pair[2],pair[1]], 12),
          limma.Rg.p = round(pv.matrix[pair[4],pair[3]], 12),
          limma.Lg.logFC = round(logFC.matrix[pair[2],pair[1]],4),
          limma.Rg.logFC = round(logFC.matrix[pair[4],pair[3]],4))
      }

    }

  }

  colnames(genepairs2) <- c("ligand_cell", "ligand_gene", "receptor_cell", "receptor_gene")


  ### parallel computation
  if (parallel.use) {
    cl <- makeCluster(cpucores)
    clusterEvalQ(cl, set.seed(2021))
    clusterExport(cl, varlist = c("samplematrix1", "samplematrix2",
                                  "mu.matrix1", "mu.matrix2",
                                  "sigma.matrixNull",
                                  "sigma.matrix1D", "sigma.matrix2D",
                                  "pv.matrix", "logFC.matrix",
                                  "threshold1", "threshold2",
                                  "const1", "const2",
                                  "adjust", "rhos",
                                  "Brep", "Brep0", "pv0",
                                  "t.test", "wilcox.test"), envir = environment())
    Rs <- cbind(genepairs2, data.frame(t(parApply(cl, genepairs, 1, comparefun))))
    stopCluster(cl)
  } else {
    set.seed(2021)
    Rs <- cbind(genepairs2, data.frame(t(apply(genepairs, 1, comparefun))))
  }

  if (!show.all) {
    Rs <- Rs[Rs$stage %in% c(1,2),]
    Rs$pvalue <- as.numeric(as.character(Rs$pvalue))
    Rs$obs.xy.diff <- as.numeric(as.character(Rs$obs.xy.diff))
    Rs$null.diff.sd <- as.numeric(as.character(Rs$null.diff.sd))
    Rs$diff.stat <- as.numeric(as.character(Rs$diff.stat))
    Rs$pvalue <- as.numeric(as.character(Rs$pvalue))
    Rs$Welch.t.stat <- as.numeric(as.character(Rs$Welch.t.stat))
    Rs$Welch.t.sd <- as.numeric(as.character(Rs$Welch.t.sd))
    Rs$Welch.t.p <- as.numeric(as.character(Rs$Welch.t.p))
    Rs$WRS.stat <- as.numeric(as.character(Rs$WRS.stat))
    Rs$WRS.p <- as.numeric(as.character(Rs$WRS.p))
    Rs$rhoNull <- as.numeric(as.character(Rs$rhoNull))
    Rs$limma.Lg.p <- as.numeric(as.character(Rs$limma.Lg.p))
    Rs$limma.Rg.p <- as.numeric(as.character(Rs$limma.Rg.p))
    Rs$limma.Lg.logFC <- as.numeric(as.character(Rs$limma.Lg.logFC))
    Rs$limma.Rg.logFC <- as.numeric(as.character(Rs$limma.Rg.logFC))
  }
  list(Rs=Rs,
       mu.matrix1=mu.matrix1, mu.matrix2=mu.matrix2,
       sigma.matrixNull=sigma.matrixNull,
       sigma.matrix1=sigma.matrix1D,
       sigma.matrix2=sigma.matrix2D,
       pv.matrix=pv.matrix,
       logFC.matrix=logFC.matrix
  )
}
