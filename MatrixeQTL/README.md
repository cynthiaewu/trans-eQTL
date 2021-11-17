# Matrix eQTL

 ## Files
 
 ### Run from step 1
 Genotype file
 
 Expression file
 
 Covariates file (optional)
 
 ## Steps
1. Run Matrix eQTL. Input genotype and expression files. Output matrix eqtl file of pairwise snp-gene pairs
   ```
   Rscript Run_Matrix_eQTL_PC_PEER_quantile_norm.r -g {genotype} -e {expression} -c {covariates} -o {matrixeqtl_file}
   ```
