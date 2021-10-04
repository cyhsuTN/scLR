# scLR
scLR: a method to test dysregulated ligand-receptor interactions

1. Download and install scLR_0.9.4.tar.gz
2. See Demo_scLR_0903.pdf for usage

Version 0.9.4: 
correct an error when setting impute.miss.celltype = 0.

Version 0.9.3:
Add parameters
- normalization. Perform normalization within each cell type (normalization = "ByCT") or for all cell types together (normalization = "ByALL").
- rhos. Estimate correlation coefficient for statistically significant lr pair (rhos = "est").
- impute.miss.celltype. Impute NA columns for missing cell types of replicates. 
