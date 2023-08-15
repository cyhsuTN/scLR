scLR: a method to test dysregulated ligand-receptor interactions
================
Chih-Yuan Hsu\
April 18, 2022

## Version

-   Version 0.9.6: updated on 04/18/2022:

1.  Correct a bug when each condition has only one sample in some LR comparison.

-   Version 0.9.5: updated on 12/30/2021:

1.  Add the nonparametric Wilcoxon rank sum test to evaluate LR pairs.
2.  Add parameters
    \_ bntest. (default = FALSE) Whether or not to perform a bivariate noraml test to check bivariate normal assumption for LR pairs.
    \_ sig.bntest. (default = 0.01) Under what significance level, the LR pairs will be tested by a bivariate noraml test.

-   Version 0.9.4:

1.  Correct an error when setting impute.miss.celltype = 0.

-   Version 0.9.3:

1.  Add parameters
    \_ normalization. Perform normalization within each cell type (normalization = "ByCT") or for all cell types together (normalization = "ByALL").
    \_ rhos. Estimate correlation coefficient for statistically significant lr pair (rhos = "est").
    \_ impute.miss.celltype. Impute NA columns for missing cell types of replicates.

## Usage of scLR

1.  Download and install scLR\_0.9.6.tar.gz


A simulated data to compare the products of the expressions of LR pairs between two conditions: TX1 and TX2. 3 samples for each condition. For each sample, there are 1000 genes and 100 cells (5 cell types, 20 cells each). 10 ligand-receptor gene pairs across 5 cell types are compared between TX1 and TX2.

Three inputs are required: countmatrix (gene names required), cellinfo, lrpairs.sample

``` r
 library(scLR)
 set.seed(2021)
 G <- 1000; n <- 600 # To create a simulated data consisting of 1000 genes and 600 cells
 # Data are generated from NB distribution
 mu1 <- rgamma(G, shape = 2, rate = 2)
 NB_cell <- function(j) rnbinom(G, size = 0.1, mu = mu1)
 countmatrix <- as(sapply(1:n, NB_cell), "sparseMatrix") 
 # 1000 pairs of LR names are sampled from lrpairs0 (built in scLR),
 # which is a list of LR commonly compared.
 genenames <- unique(unlist(lrpairs0))
 rownames(countmatrix) <- genenames[sample(1:length(genenames),1000)]
 
 # Information for all cells
 cellinfo <- data.frame(sampleID = factor(paste0("s", rep(1:6, each=100))),
                       condition = factor(paste0("tx", rep(1:2, each=300))),
                       cellcluster = factor(paste0("cc", rep(rep(1:5, each=20), 6))) )

 # Names of 10 ligand-receptor pairs which will be compared
 lrpairs.sample <- data.frame(lrpairs0[sample(1:200, 10),])
```

Formats of 3 inputs:

``` r
 countmatrix[1:6, 1:10] # format of countmatrix, gene expressions (the first 6 genes and 10 cells)
```

    ## 6 x 10 sparse Matrix of class "dgCMatrix"
    ##                             
    ## NID1    15 . . . . 2 . . 1 .
    ## IAPP     . . . 2 . . . . . .
    ## MRGPRX1  1 . 1 4 . . . . . .
    ## GPR182   . 3 . . . . . . . 3
    ## CD58     1 1 . . . . . 7 4 5
    ## SLC45A3  . 1 . . . . . 1 1 .

``` r
 head(cellinfo, 10) # format of cellinfo, categories of cells (the first 10 cells)
```

    ##    sampleID condition cellcluster
    ## 1        s1       tx1         cc1
    ## 2        s1       tx1         cc1
    ## 3        s1       tx1         cc1
    ## 4        s1       tx1         cc1
    ## 5        s1       tx1         cc1
    ## 6        s1       tx1         cc1
    ## 7        s1       tx1         cc1
    ## 8        s1       tx1         cc1
    ## 9        s1       tx1         cc1
    ## 10       s1       tx1         cc1

``` r
 head(lrpairs.sample, 10) # format of lrpairs.sample, LR gene pairs to be compared.
```

    ##     ligand receptor
    ## 200   ASIP     MC1R
    ## 94    AHSG     INSR
    ## 126   APLN   ADRA2A
    ## 56    ADM2    GPR84
    ## 173   APOE    VLDLR
    ## 23  ADAM23    ITGA5
    ## 33   ADAM9    ITGAV
    ## 123  ANXA1     FPR3
    ## 80    AGRN     LRP2
    ## 129   APLN   MTNR1A

Output:

``` r
 output <- scLR(countmatrix, cellinfo, lrpairs.sample, low.filter = 1,
                parallel.use = FALSE) # Do parallel computation if parallel.use = TRUE.
```

    ## converting counts to integer mode

``` r
 head(output$Rs[,1:11], 10)
```

    ##    lr.cell.name lr.gene.name obs.xy.diff null.diff.sd pvalue stage     adj.p
    ## 1       cc1-cc1    AHSG-INSR     -5.2920       4.9666  0.280     1 0.9147727
    ## 2       cc1-cc1   ADM2-GPR84     -2.4461       5.9463  0.680     1 0.9147727
    ## 3       cc1-cc1   APOE-VLDLR      6.8263       5.6654  0.230     1 0.9147727
    ## 4       cc1-cc1    AGRN-LRP2      1.2687       4.4509  0.800     1 0.9147727
    ## 5       cc1-cc2    AHSG-INSR     -3.6669       5.2208  0.505     1 0.9147727
    ## 6       cc1-cc2   ADM2-GPR84      2.2378       6.1187  0.655     1 0.9147727
    ## 7       cc1-cc2   APOE-VLDLR     -3.6713       4.8951  0.425     1 0.9147727
    ## 8       cc1-cc2    AGRN-LRP2      4.3512       5.3039  0.385     1 0.9147727
    ## 9       cc1-cc3    AHSG-INSR     -2.1883       5.1272  0.695     1 0.9147727
    ## 10      cc1-cc3   ADM2-GPR84      1.8088       5.6154  0.750     1 0.9147727
    ##    Welch.t.stat Welch.t.sd Welch.t.p Welch.t.adj.p
    ## 1       -0.9930     5.3294  0.391501     0.9267978
    ## 2       -0.4570     5.3522  0.677112     0.9267978
    ## 3        1.3982     4.8821  0.269773     0.9267978
    ## 4        0.3150     4.0270  0.768918     0.9267978
    ## 5       -1.8562     1.9755  0.145036     0.9267978
    ## 6        0.4669     4.7929  0.684686     0.9267978
    ## 7       -0.5773     6.3589  0.597654     0.9267978
    ## 8        0.7012     6.2055  0.527396     0.9267978
    ## 9       -0.2922     7.4899  0.784703     0.9267978
    ## 10       0.2549     7.0966  0.819315     0.9267978

## Assume sample s2 does NOT have cell type cc3

``` r
idx.remove <- which(cellinfo$sampleID=="s2" & cellinfo$cellcluster=="cc3")
cellinfo1 <- cellinfo[-idx.remove,]; cellinfo1 <- droplevels(cellinfo1)
countmatrix1 <- countmatrix[,-idx.remove]
```

Output: (The missing cell type cc3 in sample s2 will be replaced by NA)

``` r
output <- scLR(countmatrix1, cellinfo1, lrpairs.sample, low.filter = 1, 
               parallel.use = FALSE, impute.miss.celltype = NA)
```

    ## converting counts to integer mode

``` r
head(output$Rs[,1:11], 10)
```

    ##    lr.cell.name lr.gene.name obs.xy.diff null.diff.sd pvalue stage     adj.p
    ## 1       cc1-cc1    AHSG-INSR     -5.2976       4.9828  0.285     1 0.9204545
    ## 2       cc1-cc1   ADM2-GPR84     -2.4810       6.1074  0.665     1 0.9204545
    ## 3       cc1-cc1   APOE-VLDLR      6.8325       5.7146  0.240     1 0.9204545
    ## 4       cc1-cc1    AGRN-LRP2      1.2488       4.4722  0.800     1 0.9204545
    ## 5       cc1-cc2    AHSG-INSR     -3.6596       5.2088  0.490     1 0.9204545
    ## 6       cc1-cc2   ADM2-GPR84      2.2899       6.2846  0.665     1 0.9204545
    ## 7       cc1-cc2   APOE-VLDLR     -3.6278       4.9547  0.445     1 0.9204545
    ## 8       cc1-cc2    AGRN-LRP2      4.3865       5.2130  0.355     1 0.9204545
    ## 9       cc1-cc3    AHSG-INSR      2.2369       5.9956  0.735     1 0.9204545
    ## 10      cc1-cc3   ADM2-GPR84      1.7947       6.4591  0.760     1 0.9204545
    ##    Welch.t.stat Welch.t.sd Welch.t.p Welch.t.adj.p
    ## 1       -0.9850     5.3780  0.395984     0.9310549
    ## 2       -0.4588     5.4070  0.675716     0.9310549
    ## 3        1.3734     4.9750  0.277129     0.9310549
    ## 4        0.3065     4.0749  0.774943     0.9310549
    ## 5       -1.8344     1.9950  0.148310     0.8892294
    ## 6        0.4733     4.8376  0.680849     0.9310549
    ## 7       -0.5688     6.3784  0.602728     0.9310549
    ## 8        0.7033     6.2369  0.525957     0.9310549
    ## 9        0.3007     7.4385  0.785657     0.9310549
    ## 10       0.2341     7.6670  0.830869     0.9310549

## References
Qi Liu, Chih-Yuan Hsu*, Jia Li, Yu Shyr (2022). Dysregulated Ligand-receptor interactions from single cell transcriptomics. Bioinformatics, doi:10.1093/bioinformatics/btac294. *Co-first authorship.
