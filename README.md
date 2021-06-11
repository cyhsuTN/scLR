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
