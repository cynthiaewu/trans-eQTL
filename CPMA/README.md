# CPMA

 ## Files
 
 ### Run from step 1
 Genotype file
 
 Expression file
 
 Covariates file (optional)
 
 #### OR
 
 ### Run from step 2
 
 pairwise matrix-eqtl results
 
 ## To run all steps below in one script for files obtained from the simulation framework
 
 folder = path to folder with all simulated file
 scripts_folder = path to all scripts
 method = 1 for CPMA, 2 for x-QTL, 0 to run both
 iterations = number of simulations
 -q to run matrix eqtl first to get pairwise snp-gene pairs values for CPMA and/or x-QTL
 -g to account for gene correlation; only for simulations with gene correlation
 
 ``` 
 python run_cpma_xqtl_pipeline.py -f folder -p scripts_folder -m method -i iterations -q -g 
 ```
 
 ## Steps
 
1. Run Matrix eQTL. Input genotype and expression files. Output matrix eqtl file of pairwise snp-gene pairs
   ```
   Rscript ../MatrixeQTL/Run_Matrix_eQTL_PC_PEER_quantile_norm.r -g {genotype} -e {expression} -o {matrixeqtl_file}
   ```
2. Get the zscores (T-stat) in snp by gene format from Matrix eQTL output. Input matrix eqtl results from Step 1. Output zscore file
   ```
   python get_values.py -i {matrixeqtl_file} -z {zscore_file}
   ```
3. Get the pvalues from the T-stats file. Input zscore_file from Step 2. Output pvalue file
   ```
   python tstat_to_pvalue.py -i {zscore_file} -p {pvalue_file}
   ```

4. Calculate CPMA and/or x-QTL values for each snp from pvalues file. Input pvalue file from Step 3. Output cpma or xqtl file
   ```
   python calculate_cpma.py -i {pvalues_file} -o {cpma_file}
   ```
   ```
   python ../x-QTL/calculate_mixturemodel.py -i {pvalues_file} -o {xqtl_file}
   ```
### If no gene correlation needs to be adjusted for (simulation framework with identity set to True), follow step 5
5. Obtain pvalue for each snp from cpma values. Since there is no gene correlation, cpma values can be compared to the chi distribution with df=1. Input cpma_file from Step 4. Output empirical_file with pvalue for each snp
   ```
   python ../Simulator/compute_pvalue_chidist_nocorr.py  -i {cpma_file} -o {empirical_file}
   ```
      
### If gene correlation needs to be adjusted for, follow steps 6-8
6. Calculate the eigendecomposition from the gene covariance matrix observed from the zscore file. Also obtain the mean zscores for genes. Input zscore file from Step 2. Output mean zscores (mzscores), eigenvalues (evalues_file), and eigenvectors (evectors_file)
   ```
   python calculate_cov_meanzscores_edecomposition.py -i {zscore_file} -m {mzscores} -e {evalues_file} -q {evectors_file}
   ```

7. Simulate empirical cpma values and xqtl values from normal distribution with eigendecomposition values and mean zscores. Input mean zscores (mzscores), eigenvalues (evalues_file), and eigenvectors (evectors_file) from Step 5. Output sim_file with simulated empirical null cpma and sim_mix_file with simulated empirical null xqtl values. num_sim = number of simulated cpma/xqtl values. Recommended num_sim = 500000

   ```
   python ../x-QTL/simulate_cpma-mix_chunks.py -z {mzscores} -e {evalues_file} -q {evectors_file} -o {sim_file} -m {sim_mix_file} -n {num_sim}
   ```
 
8. Compare simulated cpma with observed cpma to get an empirical pvalue for each snp. Input sim_file from Step 6 and cpma_file from Step 4. Output empirical_file with empirical pvalue for each snp
   ```
   python calculate_empirical_pvalue.py -s {sim_file} -o {cpma_file} -e {empirical_file}
   ```
*** TODO: xqtl empirical pvalues and simulate empirical null from randomly permuted snp genotype of different allele freq or bin diff allele freq instead of eigendecomposition
