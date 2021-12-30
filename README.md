# scLR
scLR: a method to test dysregulated ligand-receptor interactions

1. Download and install scLR_0.9.5.tar.gz
2. See Demo_scLR_0903.pdf for usage

Version 0.9.5: updated on 12/30/2021:\
a. Add the nonparametric Wilcoxon rank sum test to evaluate LR pairs.\
b. Add parameters\
_ bntest. (default = FALSE) Whether or not to perform a bivariate noraml test to check bivariate normal assumption for LR pairs.\
_ sig.bntest. (default = 0.01) Under what significance level, the LR pairs will be tested by a bivariate noraml test.

Version 0.9.4:\
a. Correct an error when setting impute.miss.celltype = 0.

Version 0.9.3:\
a. Add parameters\
_ normalization. Perform normalization within each cell type (normalization = "ByCT") or for all cell types together (normalization = "ByALL").
_ rhos. Estimate correlation coefficient for statistically significant lr pair (rhos = "est").
_ impute.miss.celltype. Impute NA columns for missing cell types of replicates. 
