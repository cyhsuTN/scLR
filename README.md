# scLR
scLR: a method to test dysregulated ligand-receptor interactions

1. Download and install scLR_0.8.0.tar.gz
2. See Demo_scLR_0602.html for usage

---
title: "R function scLR"
author: "C.-Y. Hsu"
date: "June 02, 2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Usage of scLR

A simulated data to compare LR pairs between two conditions: TX1 and TX2. 3 replicates each condition. 
For each replicate, 1000 genes and 100 cells (5 cell types, 20 cells each).
10 ligand-receptor gene pairs across 5 cell types are compared between TX1 and TX2.


3 inputs required: countmatrix (gene names required), cellinfo, lrpairs.sample
```{r demo scLR}
 library(scLR)
 set.seed(2021)
 G <- 1000; n <- 600 # To create a simulated data consisting of 1000 genes and 600 cells
 # Data are generated from NB distribution
 NB_cell <- function(j) rnbinom(G, size = 0.1, mu = rgamma(G, shape = 2, rate = 2))
 countmatrix <- as(sapply(1:n, NB_cell), "sparseMatrix") 
 # 1000 gene names are taken from lrpairs0 [LR pairs which are commonly compared (built in scLR)]
 genenames <- unique(unlist(lrpairs0))
 rownames(countmatrix) <- genenames[sample(1:length(genenames),1000)]
 
 # Information for all cells
 cellinfo <- data.frame(sampleID = factor(paste0("s", rep(1:6, each=100))),
                       condition = factor(paste0("tx", rep(1:2, each=300))),
                       cellcluster = factor(paste0("cc", rep(rep(1:5, each=20), 6))) )

 # Names of 10 ligand-receptor pairs which will be compared
 lrpairs.sample <- data.frame(lrpairs0[sample(1:200, 10),])
```
