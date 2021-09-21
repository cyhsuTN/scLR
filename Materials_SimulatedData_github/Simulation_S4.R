library(dplyr)
library(limma)
library(scLR)
library(Matrix)
options(scipen = 999)
setwd("C:/Users/farca/Desktop/Project_SingleCellSpatial/Data")



mu.matrix1 <- mu.matrix2 <- ave.mu.matrix <- as.matrix(read.csv("mu.matrix.simulation.csv", fileEncoding="UTF-8-BOM"))
sigma.matrix1 <- sigma.matrix2 <- as.matrix(read.csv("sd.matrix.simulation.csv", fileEncoding="UTF-8-BOM"))

#set.seed(2021)
simulation.fun <- function(signal, pv.cri=0.01, ave.mu.matrix, mu.matrix1, mu.matrix2,
                           sigma.matrix1, sigma.matrix2, zero.impute, rhos, sizes) {

  if (sizes == 5) {
    #set.seed(2021+k)
    s1 <- matrix(rnorm(length(mu.matrix1), mean=c(ave.mu.matrix), sd=c(sigma.matrix1)), ncol=11)
    s3 <- matrix(rnorm(length(mu.matrix1), mean=c(ave.mu.matrix), sd=c(sigma.matrix1)), ncol=11)
    s4 <- matrix(rnorm(length(mu.matrix1), mean=c(ave.mu.matrix), sd=c(sigma.matrix1)), ncol=11)
    s7 <- matrix(rnorm(length(mu.matrix1), mean=c(ave.mu.matrix), sd=c(sigma.matrix1)), ncol=11)
    s8 <- matrix(rnorm(length(mu.matrix1), mean=c(ave.mu.matrix), sd=c(sigma.matrix1)), ncol=11)

    s2 <- matrix(rnorm(length(mu.matrix2), mean=c(ave.mu.matrix), sd=c(sigma.matrix2)), ncol=11)
    s5 <- matrix(rnorm(length(mu.matrix2), mean=c(ave.mu.matrix), sd=c(sigma.matrix2)), ncol=11)
    s6 <- matrix(rnorm(length(mu.matrix2), mean=c(ave.mu.matrix), sd=c(sigma.matrix2)), ncol=11)
    s9 <- matrix(rnorm(length(mu.matrix2), mean=c(ave.mu.matrix), sd=c(sigma.matrix2)), ncol=11)
    s10 <- matrix(rnorm(length(mu.matrix2), mean=c(ave.mu.matrix), sd=c(sigma.matrix2)), ncol=11)


    #set.seed(2021)
    n.pairs <- 1000
    genepairs1 <- matrix(NA, nrow=n.pairs, ncol=2)
    colnames(genepairs1) <- c("cell", "gene")
    genepairs1[,c(1)] <- sample(1:11, n.pairs, replace = TRUE)
    genepairs1[,c(2)] <- sample(1:16672, n.pairs, replace = TRUE)

    DE.gene <- data.frame(genepairs1[,c(2,1)])
    DE.gene.unique <- distinct(DE.gene, gene, cell, .keep_all= TRUE)

    idx.DE <- as.matrix(DE.gene.unique)
    s2[idx.DE[1:500,]] <- s2[idx.DE[1:500,]] * signal
    s5[idx.DE[1:500,]] <- s5[idx.DE[1:500,]] * signal
    s6[idx.DE[1:500,]] <- s6[idx.DE[1:500,]] * signal
    s9[idx.DE[1:500,]] <- s9[idx.DE[1:500,]] * signal
    s10[idx.DE[1:500,]] <- s10[idx.DE[1:500,]] * signal
    s2[idx.DE[501:nrow(DE.gene.unique),]] <- s2[idx.DE[501:nrow(DE.gene.unique),]] / signal
    s5[idx.DE[501:nrow(DE.gene.unique),]] <- s5[idx.DE[501:nrow(DE.gene.unique),]] / signal
    s6[idx.DE[501:nrow(DE.gene.unique),]] <- s6[idx.DE[501:nrow(DE.gene.unique),]] / signal
    s9[idx.DE[501:nrow(DE.gene.unique),]] <- s9[idx.DE[501:nrow(DE.gene.unique),]] / signal
    s10[idx.DE[501:nrow(DE.gene.unique),]] <- s10[idx.DE[501:nrow(DE.gene.unique),]] / signal

    samplematrix1 <- cbind(s1, s3, s4, s7, s8)
    samplematrix2 <- cbind(s2, s5, s6, s9, s10)
  } else {
    #set.seed(2021+k)
    s4 <- matrix(rnorm(length(mu.matrix1), mean=c(ave.mu.matrix), sd=c(sigma.matrix1)), ncol=11)
    s7 <- matrix(rnorm(length(mu.matrix1), mean=c(ave.mu.matrix), sd=c(sigma.matrix1)), ncol=11)
    s8 <- matrix(rnorm(length(mu.matrix1), mean=c(ave.mu.matrix), sd=c(sigma.matrix1)), ncol=11)

    s6 <- matrix(rnorm(length(mu.matrix2), mean=c(ave.mu.matrix), sd=c(sigma.matrix2)), ncol=11)
    s9 <- matrix(rnorm(length(mu.matrix2), mean=c(ave.mu.matrix), sd=c(sigma.matrix2)), ncol=11)
    s10 <- matrix(rnorm(length(mu.matrix2), mean=c(ave.mu.matrix), sd=c(sigma.matrix2)), ncol=11)


    #set.seed(2021)
    n.pairs <- 1000
    genepairs1 <- matrix(NA, nrow=n.pairs, ncol=2)
    colnames(genepairs1) <- c("cell", "gene")
    genepairs1[,c(1)] <- sample(1:11, n.pairs, replace = TRUE)
    genepairs1[,c(2)] <- sample(1:16672, n.pairs, replace = TRUE)

    DE.gene <- data.frame(genepairs1[,c(2,1)])
    DE.gene.unique <- distinct(DE.gene, gene, cell, .keep_all= TRUE)

    idx.DE <- as.matrix(DE.gene.unique)
    s6[idx.DE[1:500,]] <- s6[idx.DE[1:500,]] * signal
    s9[idx.DE[1:500,]] <- s9[idx.DE[1:500,]] * signal
    s10[idx.DE[1:500,]] <- s10[idx.DE[1:500,]] * signal
    s6[idx.DE[501:nrow(DE.gene.unique),]] <- s6[idx.DE[501:nrow(DE.gene.unique),]] / signal
    s9[idx.DE[501:nrow(DE.gene.unique),]] <- s9[idx.DE[501:nrow(DE.gene.unique),]] / signal
    s10[idx.DE[501:nrow(DE.gene.unique),]] <- s10[idx.DE[501:nrow(DE.gene.unique),]] / signal

    samplematrix1 <- cbind(s4, s7, s8)
    samplematrix2 <- cbind(s6, s9, s10)
  }


  genepairs2 <- matrix(NA, nrow=4*n.pairs, ncol=2)
  colnames(genepairs2) <- c("cell", "gene")
  genepairs2[,c(1)] <- sample(1:11, 4*n.pairs, replace = TRUE)
  genepairs2[,c(2)] <- sample(1:16672, 4*n.pairs, replace = TRUE)
  non.DE.gene <- data.frame(genepairs2[,c(2,1)])

  geneall <- rbind(DE.gene.unique, non.DE.gene)
  geneall.unique <- distinct(geneall, gene, cell, .keep_all= TRUE)

  l2   <- as.matrix(geneall.unique[sample(1:500, 500, replace = T),])
  r2   <- as.matrix(geneall.unique[sample(501:nrow(DE.gene.unique), 500, replace = T),])

  idx.l3 <- sample((nrow(DE.gene.unique)+1):nrow(geneall.unique), 1500, replace = T)
  idx.r3 <- sample(setdiff((nrow(DE.gene.unique)+1):nrow(geneall.unique), idx.l3), 1500, replace = T)
  l3   <- as.matrix(geneall.unique[idx.l3,])
  r3   <- as.matrix(geneall.unique[idx.r3,])

  lr2 <- distinct(data.frame(cbind(l2, r2)), gene, cell, gene.1, cell.1, .keep_all= TRUE); dim(lr2)
  lr3 <- distinct(data.frame(cbind(l3, r3)), gene, cell, gene.1, cell.1, .keep_all= TRUE); dim(lr3)

  genepairs <- (rbind(lr2, lr3))[,c(2,1,4,3)]
  rownames(genepairs) <- 1:nrow(genepairs)

  m.DE <- nrow(lr2)


  samplematrix1[samplematrix1<0] <- 0
  samplematrix2[samplematrix2<0] <- 0

  rownames(samplematrix1) <- rownames(samplematrix2) <- paste0("g",1:nrow(samplematrix1))
  colnames(samplematrix1) <- colnames(samplematrix2) <- paste0("c",rep(1:11, ncol(samplematrix1)/11))


  ###
  ### Limma
  ###
  N1 <- ncol(samplematrix1)/11
  N2 <- ncol(samplematrix2)/11
  const1 <- (0:(N1 - 1)) * 11
  const2 <- (0:(N2 - 1)) * 11
  bb <- NULL
  for (k in 1:11) {
    idx.col1 <- k + const1
    idx.col2 <- k + const2

    yy <- cbind(samplematrix1[, idx.col1], samplematrix2[, idx.col2])
    designmatrix <- cbind(Grp1=1, Grp2vs1=c(rep(0,N1),rep(1,N2)))

    fit <- lmFit(yy, designmatrix)
    fit <- eBayes(fit)
    aa <- topTable(fit, coef=2, p.value=pv.cri, num=1000)
    if (nrow(aa)>0) {
      bb <- rbind(bb, cbind(celltype=k, gene=pmatch(rownames(aa), rownames(samplematrix1)), aa))
    }
  }
  bbgc <- paste0(bb$gene,"-",bb$celltype)
  lr2.gc12 <- paste0(lr2[,1],"-",lr2[,2])
  lr2.gc34 <- paste0(lr2[,3],"-",lr2[,4])
  lr3.gc12 <- paste0(lr3[,1],"-",lr3[,2])
  lr3.gc34 <- paste0(lr3[,3],"-",lr3[,4])


  idx.limma.sig.gene.pairs <- unique( c(which(lr2.gc12 %in% bbgc | lr2.gc34 %in% bbgc),
                                        which(lr3.gc12 %in% bbgc | lr3.gc34 %in% bbgc) + nrow(lr2)) )


  #start_time <- Sys.time()
  #set.seed(0328)
  RR1.all <- scXY2(samplematrix1, samplematrix2, genepairs, No_celltypes=11, Brep=1000,
                  Brep0=200, pv0=0.05, p.filter = 0.10, show.all = F, parallel.use = T, cpucores = 2, adjust = 3,
                  zero.impute = zero.impute, rhos = rhos)

  #end_time <- Sys.time()
  #end_time - start_time

  RR <- RR1.all$Rs
  RR$adj.p <- round(p.adjust(as.numeric(as.character(RR$pvalue)), method="BH"), 12)
  RR$adj.t.p <- round(p.adjust(as.numeric(as.character(RR$Welch.t.p)), method="BH"), 12)


  sig.pair.name.scXY <- rownames(RR)[which(RR$adj.p < pv.cri)]
  num.p.scXY <- length(sig.pair.name.scXY)
  num.true.p.scXY <- 0
  precesion.scXY <- 0
  recall.scXY <- 0
  FalseP.scXY <- length(sig.pair.name.scXY)/nrow(genepairs)
  F1.scXY <- NA
  one.p.scXY <- NA

  sig.pair.name.t <- rownames(RR)[which(RR$adj.t.p < pv.cri)]
  num.p.t <- length(sig.pair.name.t)
  num.true.p.t <- 0
  precesion.t <- 0
  recall.t <- 0
  FalseP.t <- length(sig.pair.name.t)/nrow(genepairs)
  F1.t <- NA
  one.p.t <- NA

  sig.pair.name.limma <- rownames(RR)[idx.limma.sig.gene.pairs]
  num.p.limma <- length(sig.pair.name.limma)
  num.true.p.limma <- 0
  precesion.limma <- 0
  recall.limma <- 0
  FalseP.limma <- length(sig.pair.name.limma)/nrow(genepairs)
  F1.limma <- NA
  one.p.limma <- NA

  ROC.scXY <- function(pvcri) {
    sig.pair.name.scXY <- rownames(RR)[which(RR$adj.p < pvcri)]
    recall.scXY <- 0
    FalseP.scXY <- length(sig.pair.name.scXY)/nrow(genepairs)

    num.p.scXY <- length(sig.pair.name.scXY)
    num.true.p.scXY <- 0
    precesion.scXY <- 0
    F1.scXY <- NA
    c(TP=recall.scXY, FP=FalseP.scXY, F1=F1.scXY)
  }

  ROC.t <- function(pvcri) {
    sig.pair.name.t <- rownames(RR)[which(RR$adj.t.p < pvcri)]
    recall.t <- 0
    FalseP.t <- length(sig.pair.name.t)/nrow(genepairs)

    num.p.t <- length(sig.pair.name.t)
    num.true.p.t <- 0
    precesion.t <- 0
    F1.t <- NA
    c(TP=recall.t, FP=FalseP.t, F1=F1.t)
  }

  ROC.limma <- function(pvcri) {
    N1 <- ncol(samplematrix1)/11
    N2 <- ncol(samplematrix2)/11
    const1 <- (0:(N1 - 1)) * 11
    const2 <- (0:(N2 - 1)) * 11
    bb <- NULL
    for (k in 1:11) {
      #k <- 3
      idx.col1 <- k + const1
      idx.col2 <- k + const2

      yy <- cbind(samplematrix1[, idx.col1], samplematrix2[, idx.col2])
      designmatrix <- cbind(Grp1=1, Grp2vs1=c(rep(0,N1),rep(1,N2)))

      fit <- lmFit(yy, designmatrix)
      fit <- eBayes(fit)
      aa <- topTable(fit, coef=2, p.value=pvcri, num=1000)
      if (nrow(aa)>0) {
        bb <- rbind(bb, cbind(celltype=k, gene=pmatch(rownames(aa), rownames(samplematrix1)), aa))
      }
    }
    bbgc <- paste0(bb$gene,"-",bb$celltype)
    lr2.gc12 <- paste0(lr2[,1],"-",lr2[,2])
    lr2.gc34 <- paste0(lr2[,3],"-",lr2[,4])
    lr3.gc12 <- paste0(lr3[,1],"-",lr3[,2])
    lr3.gc34 <- paste0(lr3[,3],"-",lr3[,4])


    idx.limma.sig.gene.pairs <- unique( c(which(lr2.gc12 %in% bbgc | lr2.gc34 %in% bbgc),
                                          which(lr3.gc12 %in% bbgc | lr3.gc34 %in% bbgc) + nrow(lr2)) )

    sig.pair.name.limma <- rownames(RR)[idx.limma.sig.gene.pairs]
    recall.limma <- 0
    FalseP.limma <- length(sig.pair.name.limma)/nrow(genepairs)

    num.p.limma <- length(sig.pair.name.limma)
    num.true.p.limma <- 0
    precesion.limma <- 0
    F1.limma <- NA
    c(TP=recall.limma, FP=FalseP.limma, F1=F1.limma)
  }

  ForROC <- t(rbind(sapply(c(0.001, 0.005, 0.01, 0.025, 0.05, 0.1), # seq(0.12,0.98,0.02), 0.999),
                           function(x) ROC.scXY(x)),
                    sapply(c(0.001, 0.005, 0.01, 0.025, 0.05, 0.1), # seq(0.12,0.98,0.02), 0.999),
                           function(x) ROC.t(x)),
                    sapply(c(0.001, 0.005, 0.01, 0.025, 0.05, 0.1), # seq(0.12,0.98,0.02), 0.999),
                           function(x) ROC.limma(x)) ))

  list(Result = c(c(num.p.scXY=num.p.scXY, num.true.p.scXY=num.true.p.scXY, precesion.scXY=precesion.scXY,
                    recall.scXY=recall.scXY, FalseP.scXY=FalseP.scXY, F1.scXY=F1.scXY, one.p.scXY=one.p.scXY),
                  c(num.p.t=num.p.t, num.true.p.t=num.true.p.t, precesion.t=precesion.t,
                    recall.t=recall.t, FalseP.t=FalseP.t, F1.t=F1.t, one.p.t=one.p.t),
                  c(num.p.limma=num.p.limma, num.true.p.limma=num.true.p.limma, precesion.limma=precesion.limma,
                    recall.limma=recall.limma, FalseP.limma=FalseP.limma, F1.limma=F1.limma, one.p.limma=one.p.limma)),
       ForROC = ForROC)
}


zero.impute <- TRUE
#signal <- 2; sizes <- 5; seed <- 405 # for 5 samples
signal <- 2; sizes <- 3; seed <- 505 # for 3 samples
rhos <- 0
sim.results <- sapply(1:10, function(i) {
  set.seed(seed+i)
  simulation.fun(signal, pv.cri=0.01, ave.mu.matrix, mu.matrix1, mu.matrix2,
                 sigma.matrix1, sigma.matrix2, zero.impute, rhos, sizes)
})

(sim.results0504 <- sapply(sim.results[1,], function(x) x))
write.csv(sim.results0504, file = 'sim.results_3v3_1_S4.csv', row.names = T)
#write.csv(sim.results0504, file = 'sim.results_5v5_1_S4.csv', row.names = T)

(ROC.ave <- matrix(rowMeans( sapply(sim.results[2,], function(x) x), na.rm=T), ncol = 9))
colnames(ROC.ave) <- c("TP.scXY", "FP.scXY", "F1.scXY",
                       "TP.t", "FP.t", "F1.t",
                       "TP.limma", "FP.limma", "F1.limma")
write.csv(ROC.ave, file = 'sim.results_3v3_2_S4.csv', row.names = T)
#write.csv(ROC.ave, file = 'sim.results_5v5_2_S4.csv', row.names = T)





