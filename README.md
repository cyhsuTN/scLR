# scLR
scLR: a method to test dysregulated ligand-receptor interactions

1. Download and install scLR_0.9.2.tar.gz
2. See Demo_scLR_0903.pdf for usage

Version 0.9.2:
Add parameters
- normalization. Perform normalization within each cell type (normalization = "ByCT") or for all cell types together (normalization = "ByALL").
- rhos. Estimate correlation coefficient for each lr pair (rhos = "est").
- impute.miss.celltype. Impute NA columns for missing cell types of replicates. 
