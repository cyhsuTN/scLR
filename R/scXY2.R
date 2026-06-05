#' Test of XY at very small sample sizes (sumExp matrix)
#' @description A function to test difference of XY between two samples with very small sizes (sumExp matrix).
#' @importFrom limma lmFit eBayes
#' @importFrom stats rchisq quantile runif t.test sd bw.SJ pnorm cor.test wilcox.test
#' @importFrom Matrix head tail
#' @importFrom parallel makeCluster clusterExport parApply stopCluster clusterEvalQ


#' @param samplematrix1 Input. Data matrix/dataframe (normal-distribution-like) for condition 1 (genes by replicates * number of cell types).
#' @param samplematrix2 Input. Data matrix/dataframe (normal-distribution-like) for condition 2 (genes by replicates * number of cell types).
#' @param genepairs  Input. A K by 4 matrix/dataframe to indicate which K ligand-receptor pairs are compared between two conditions
#' It may be inputted by either index or name (See the example below for details of input).
#' @param No_celltypes  Input. Number of cell types.
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
#' @param adjust Input (optional). To adjust estimated bandwidth h (bw.SJ * adjust).
#' The default of \emph{adjust} is 3.
#' @param rhos Input (optional). The correlation coefficients in ligand-receptor pairs.
#' The default is 0, assuming zero correlations for all ligand-receptor pairs.
#' Users can assume an overall correlation value (between 0 and 1) for all ligand-receptor pairs, e.g, rhos = 0.1;
#' or set rhos = "est", the Pearson correlation coefficient estimates are used for the significant correlations
#' (p-values < 0.01 determined by R function cor.test), and 0 for the non-significant correlations.



#' @return \item{obs.xy.diff}{Observed difference of mean of XY between two conditions (condition 2 minus condition 1)}
#' @return \item{null.diff.sd}{Standard deviation of mean difference of XY under null hypothesis}
#' @return \item{diff.stat}{obs.xy.diff is divided by null.diff.sd}
#' @return \item{pvalue}{P-value}
#' @return \item{stage}{Which stage gene pairs pass to}
#' @return \item{Welch.t.*}{Statistic, pooled standard error, and p-value of Welch's t-test}
#' @return \item{WRS.*}{Statistic and p-value of Wilcoxon Rank Sum test}
#' @return \item{mu.matrix1}{Sample mean for each gene in condition 1}
#' @return \item{mu.matrix2}{Sample mean for each gene in condition 2}
#' @return \item{sigma.matrix1}{Posterior residual standard deviation for each gene in condition 1}
#' @return \item{sigma.matrix2}{Posterior residual standard deviation for each gene in condition 2}


#' @examples # An example to compare 3 and 5 replicates in TX1 and TX2, respectively.
#' @examples # 100 genes and 10 cell types for each replicate.
#' @examples # 5 ligand-receptor gene pairs are compared between replicates in TX1 and TX2.
#'
#' @examples set.seed(2025)
#' @examples ## Matrix format inputed must be as follows:
#' @examples samplematrix1 <- cbind(matrix(rnorm(100*10, mean=1.0,sd=0.3),100,10), # replicate 1 (10 cell types)
#' @examples                        matrix(rnorm(100*10, mean=1.0,sd=0.3),100,10), # replicate 2 (10 cell types)
#' @examples                        matrix(rnorm(100*10, mean=1.0,sd=0.3),100,10)) # replicate 3 (10 cell types)
#' @examples samplematrix2 <- cbind(matrix(rnorm(100*10, mean=0.65,sd=0.3),100,10), # replicate 1 (10 cell types)
#' @examples                        matrix(rnorm(100*10, mean=0.65,sd=0.3),100,10), # replicate 2 (10 cell types)
#' @examples                        matrix(rnorm(100*10, mean=0.65,sd=0.3),100,10), # replicate 3 (10 cell types)
#' @examples                        matrix(rnorm(100*10, mean=0.65,sd=0.3),100,10), # replicate 4 (10 cell types)
#' @examples                        matrix(rnorm(100*10, mean=0.65,sd=0.3),100,10)) # replicate 5 (10 cell types)
#' @examples No_celltypes <- 10 ## number of cell types
#'
#' @examples colnames(samplematrix1) <- rep(paste0("cell",1:No_celltypes), 3)
#' @examples colnames(samplematrix2) <- rep(paste0("cell",1:No_celltypes), 5)
#' @examples rownames(samplematrix1) <- rownames(samplematrix2) <- paste0("gene",1:nrow(samplematrix1))
#'
#' @examples ## Set low-expressed genes for some cell types of replicates
#' @examples samplematrix1[82, 4+(0:(3 - 1)) * No_celltypes] <- c(0, 0, 0)
#' @examples #samplematrix2[82, 4+(0:(5 - 1)) * No_celltypes] <- c(0.01, 0, 0.02, 0.1, 0)
#'
#' @examples ## Set NA (not detected) for some cell types of replicates
#' @examples samplematrix1[,c(4,12,27,28)] <- NA
#' @examples samplematrix2[,c(1,13,24,25,33,36,39,47)] <- NA
#'
#' @examples ## lignad-receptor gene pair index
#' @examples genepairs <- matrix(NA, nrow=5, ncol=4)
#' @examples colnames(genepairs) <- c("ligand_cell", "ligand_gene",
#' @examples                          "receptor_cell", "receptor_gene")
#' @examples genepairs[,c(1,3)] <- sample(1:10, 2*5)
#' @examples genepairs[,c(2,4)] <- sample(1:100, 2*5)
#' @examples ## output
#' @examples output1 <- scXY2(samplematrix1, samplematrix2, genepairs, No_celltypes, Brep=1000)
#' @examples output1$Rs
#'
#' @examples ## Or use ligand-receptor gene pair name
#' @examples lcname <- c("cell8", "cell4", "cell1", "cell5", "cell10")
#' @examples lgname <- c("gene85", "gene82", "gene7", "gene94", "gene58")
#' @examples rcname <- c("cell6", "cell3", "cell2", "cell7", "cell9")
#' @examples rgname <- c("gene22", "gene57", "gene26", "gene77", "gene14")
#' @examples pairname <- data.frame(cbind(lcname, lgname, rcname, rgname))
#' @examples ## output
#' @examples output2 <- scXY2(samplematrix1, samplematrix2, pairname, No_celltypes, Brep=1000)
#' @examples output2$Rs

#' @export
#'


scXY2 <- function(samplematrix1, samplematrix2, genepairs, No_celltypes,
                  Brep = 1000, Brep0 = 200, pv0 = 0.05, p.filter = 0.1, show.all = FALSE,
                  parallel.use = FALSE, cpucores = 2, zero.impute = TRUE, adjust = 3, rhos = 0) {

  ### convert to matrix
  if(!is.matrix(samplematrix1)) {samplematrix1 <- as.matrix(samplematrix1)}
  if(!is.matrix(samplematrix2)) {samplematrix2 <- as.matrix(samplematrix2)}
  if(!is.matrix(genepairs)) {genepairs <- as.matrix(genepairs)}

  ### logic check
  if((ncol(samplematrix1) %% No_celltypes)!=0 | (ncol(samplematrix2) %% No_celltypes)!=0) {
    stop("The number of columns of datamatrix is not a multiplier of No_celltypes")
  } else {
    N1 <- ncol(samplematrix1)/No_celltypes
    N2 <- ncol(samplematrix2)/No_celltypes
    if (N1 <= 1 & N2 <= 1) {
      stop("The numbers of replicates/subjects for each condition must be larger than 1")
    }
  }

  ### pairname convert to index
  genepairs2 <- genepairs

  if (!is.numeric(genepairs)) {

    cellname <- unique(colnames(samplematrix1))
    genename <- rownames(samplematrix1)
    if(!all(genepairs[,c(1,3)] %in% cellname)) {
      stop("Some cell names in pairname cannot be found in samplematrix")
    }
    if(!all(genepairs[,c(2,4)] %in% genename)) {
      stop("Some gene names in pairname cannot be found in samplematrix")
    }

    genepairs <-   cbind(pmatch(c(genepairs[,1]), cellname, duplicates.ok = T),
                         pmatch(c(genepairs[,2]), genename, duplicates.ok = T),
                         pmatch(c(genepairs[,3]), cellname, duplicates.ok = T),
                         pmatch(c(genepairs[,4]), genename, duplicates.ok = T))
  }


  mu.matrix1 <- sigma.matrixNull <- sigma.matrix1D <- matrix(NA, nrow = nrow(samplematrix1), ncol = No_celltypes)
  mu.matrix2 <- sigma.matrix2D <- matrix(NA, nrow = nrow(samplematrix2), ncol = No_celltypes)
  #pv.matrix <- logFC.matrix <- matrix(NA, nrow = nrow(samplematrix1), ncol = No_celltypes)
  const1 <- (0:(N1 - 1)) * No_celltypes
  const2 <- (0:(N2 - 1)) * No_celltypes

  #min1 <- min(samplematrix1[samplematrix1>0], na.rm = T)/2
  #min2 <- min(samplematrix2[samplematrix2>0], na.rm = T)/2
  min1 <- apply(samplematrix1, 2, function(x) { ifelse(all(is.na(x)), NA, min(x[x>0], na.rm = T)/4) })
  min2 <- apply(samplematrix2, 2, function(x) { ifelse(all(is.na(x)), NA, min(x[x>0], na.rm = T)/4) })

  ### run limma to estimate mean and posterior residual standard deviations
  set.seed(2021)
  for (k in 1:No_celltypes) {
    #k <- 1
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
      #pv.matrix[, k]    <-  as.numeric(fit0$p.value[,2])
    } else {
      sigma.matrixNull[, k] <- NA
      #pv.matrix[, k]    <-  NA
    }

    mu.matrix1[, k] <- as.numeric(fit0$coefficients[,1])
    mu.matrix2[, k] <- as.numeric(rowSums(fit0$coefficients))
    #logFC.matrix[, k] <-  as.numeric(fit0$coefficients[,2])

    ######
    fit1 <- lmFit(samplematrix1[, idx.col1])
    fit2 <- lmFit(samplematrix2[, idx.col2])

    fit1$sigma[which(fit1$sigma==0)] <- 1E-8
    fit2$sigma[which(fit2$sigma==0)] <- 1E-8

    if (fit1$df.residual[1]!=0) {
      fit1 <- eBayes(fit1)
      sigma.matrix1D[, k] <- sqrt(fit1$s2.post)
    } else {
      sigma.matrix1D[, k] <- 0
    }

    if (fit2$df.residual[1]!=0) {
      fit2 <- eBayes(fit2)
      sigma.matrix2D[, k] <- sqrt(fit2$s2.post)
    } else {
      sigma.matrix2D[, k] <- 0
    }

  }

  ### Obtain the thresholds for each cell type in two conditions (after normalization)
  threshold1 <- apply(mu.matrix1, 2, quantile, probs = p.filter)
  threshold2 <- apply(mu.matrix2, 2, quantile, probs = p.filter)

  comparefun <- function(pair) {
    #pair <- genepairs[1933,]
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
        Welch.t.stat = "-", Welch.t.sd = "-", Welch.t.p = "-", WRS.stat = "-", WRS.p = "-")
        #limma.Lg.p = "-", limma.Rg.p = "-", limma.Lg.logFC = "-", limma.Rg.logFC = "-")
    } else if (any(is.na(c(sigma.matrixNull[pair[2], pair[1]], sigma.matrixNull[pair[4], pair[3]])))) {
      c(obs.xy.diff = "NA sigma", null.diff.sd = "-", diff.stat = "-",
        pvalue = "-", stage = "-", rhoNull = "-",
        Welch.t.stat = "-", Welch.t.sd = "-", Welch.t.p = "-", WRS.stat = "-", WRS.p = "-")
        #limma.Lg.p = "-", limma.Rg.p = "-", limma.Lg.logFC = "-", limma.Rg.logFC = "-")
    } else if (cond1 | cond2) {
      c(obs.xy.diff = paste0("Lowly expressed"), null.diff.sd = "-", diff.stat = "-",
        pvalue = "-", stage = "-", rhoNull = "-",
        Welch.t.stat = "-", Welch.t.sd = "-", Welch.t.p = "-", WRS.stat = "-", WRS.p = "-")
      #limma.Lg.p = "-", limma.Rg.p = "-", limma.Lg.logFC = "-", limma.Rg.logFC = "-")
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
      #rho1 <- 0
      #rho2 <- 0

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
        #pvalue <-  mean(abs(c(Ts0,Ts)) > abs(obs.xy.diff), na.rm=T)
        null.diff.sd <- sd(c(Ts0,Ts))
        diff.stat <- obs.xy.diff/null.diff.sd
        stage <- 2
        c(obs.xy.diff = round(obs.xy.diff,4), null.diff.sd = round(null.diff.sd,4), diff.stat = round(diff.stat,4),
          pvalue = round(pvalue,12), stage = stage, rhoNull = rho.bar,
          Welch.t.stat = Welch.t$Welch.t.stat, Welch.t.sd = Welch.t$Welch.t.sd, Welch.t.p = Welch.t$Welch.t.p,
          WRS.stat = WRS.test$WRS.stat, WRS.p = WRS.test$WRS.p)
          #limma.Lg.p = round(pv.matrix[pair[2],pair[1]], 12),
          #limma.Rg.p = round(pv.matrix[pair[4],pair[3]], 12),
          #limma.Lg.logFC = round(logFC.matrix[pair[2],pair[1]],4),
          #limma.Rg.logFC = round(logFC.matrix[pair[4],pair[3]],4))
      } else {
        null.diff.sd <- sd(Ts0)
        diff.stat <- obs.xy.diff/null.diff.sd
        stage <- 1
        c(obs.xy.diff = round(obs.xy.diff,4), null.diff.sd = round(null.diff.sd,4), diff.stat = round(diff.stat,4),
          pvalue = pvalue0, stage = stage, rhoNull = rho.bar,
          Welch.t.stat = Welch.t$Welch.t.stat, Welch.t.sd = Welch.t$Welch.t.sd, Welch.t.p = Welch.t$Welch.t.p,
          WRS.stat = WRS.test$WRS.stat, WRS.p = WRS.test$WRS.p)
          #limma.Lg.p = round(pv.matrix[pair[2],pair[1]], 12),
          #limma.Rg.p = round(pv.matrix[pair[4],pair[3]], 12),
          #limma.Lg.logFC = round(logFC.matrix[pair[2],pair[1]],4),
          #limma.Rg.logFC = round(logFC.matrix[pair[4],pair[3]],4))
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
                                  "sigma.matrix1", "sigma.matrix2",
                                  #"pv.matrix", "logFC.matrix",
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
    #Rs$limma.Lg.p <- as.numeric(as.character(Rs$limma.Lg.p))
    #Rs$limma.Rg.p <- as.numeric(as.character(Rs$limma.Rg.p))
    #Rs$limma.Lg.logFC <- as.numeric(as.character(Rs$limma.Lg.logFC))
    #Rs$limma.Rg.logFC <- as.numeric(as.character(Rs$limma.Rg.logFC))
  }
  list(Rs=Rs,
       mu.matrix1=mu.matrix1, mu.matrix2=mu.matrix2,
       sigma.matrixNull=sigma.matrixNull,
       sigma.matrix1=sigma.matrix1D,
       sigma.matrix2=sigma.matrix2D
       #pv.matrix=pv.matrix,
       #logFC.matrix=logFC.matrix
       )
}

