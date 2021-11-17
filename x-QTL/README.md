# x-QTL

 ## Files
 
 ### Run from step 1
 Genotype file
 
 Expression file
 
 Covariates file (optional)
 
 #### OR
 
 ### Run from step 2
 
 pairwise matrix-eqtl results
 
 ## Steps
 
1. Run Matrix eQTL. Input genotype and expression files. Output matrix eqtl file of pairwise snp-gene pairs
   ```
   Rscript ../MatrixeQTL/Run_Matrix_eQTL_PC_PEER_quantile_norm.r -g {genotype} -e {expression} -o {matrixeqtl_file}
   ```
2. Get the zscores (T-stat) in snp by gene format from Matrix eQTL output. Input matrix eqtl results from Step 1. Output zscore file
   ```
   python ../CPMA/get_values.py -i {matrixeqtl_file} -z {zscore_file}
   ```
3. Get the pvalues from the T-stats file. Input zscore_file from Step 2. Output pvalue file
   ```
   python ../CPMA/tstat_to_pvalue.py -i {zscore_file} -p {pvalue_file}
   ```

4. Calculate x-QTL values for each snp from pvalues file. Input pvalue file from Step 3. Output xqtl file with x-QTL values and predicted T for each snp
   ```
   python ../x-QTL/calculate_mixturemodel.py -i {pvalues_file} -o {xqtl_file}
