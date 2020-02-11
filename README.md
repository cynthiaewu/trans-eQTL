# trans-eQTL

1. Preprocess the genotype and expression file to get intersecting samples
   - preprocess_genotypefile.py -g input_genotype_file -e input_expression_file -c input_covariates_file -o out_genotype_file -p out_expression_file -q out_covariates_file
2. Filter Genotype file to get only snps in coding regions
   - getCodingRegions.py -i input_preprocess_genotype_file -o output_file -c chr_number
3. Run Matrix eQTL
   - gene-SNP_pairs.R input_coding_genotype_file input_intersect_expression_file pca_file output
4. Get the zscores (T-stat) and pvalues from Matrix eQTL output
   - get_pvalues.py
5. Calculate CPMA values for each snp from pvalues file
   - calculate_cpma.py -i input_pvalues_file -o out_cpma_file
  
    To check the distribution of cpma values follows the chi distribution of df=1, we generated random pvalues from a normal distribution and calculated cpma values for these
   - generate_random_pvalue.py 
6. Calculate the gene covariance matrix with zscores file
   - get_cov_matrix.py 
7. Calculate the eigendecomposition and mean zscores
   - calculate_edecomposition_mzscores.py -z input_zscores_file -c input_cov_matrix_file -m mean_zscores_out_file -e eigenvalues_out_file -q eigenvectors_out_file
8. Simulate zscores from normal distribution with gene covariance matrix.
   - simulate_zscores.py -z input_mean_zscores_file -e input_eigenvalues_file -q input_eigenvectors_file -o out_sim_file -n num_simulations
   - get_sim1000_zscores.py 

