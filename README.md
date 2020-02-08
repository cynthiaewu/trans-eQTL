# trans-eQTL

1. Filter Genotype file to get only snps in coding regions
   - getCodingRegions.py -i input_genotype_file -o output_file -c chr_number
2. Run Matrix eQTL
   - gene-SNP_pairs.R 
3. Get the zscores (T-stat) and pvalues from Matrix eQTL output
   - get_pvalues.py
4. Calculate CPMA values for each snp from pvalues file
   - calculate_cpma.py -i input_pvalues_file -o out_cpma_file
  
    To check the distribution of cpma values follows the chi distribution of df=1, we generated random pvalues from a normal distribution and calculated cpma values for these
   - generate_random_pvalue.py 
5. Calculate the gene covariance matrix with zscores file
   - get_cov_matrix.py 
6. Calculate the eigendecomposition and mean zscores
   - calculate_edecomposition_mzscores.py -z input_zscores_file -c input_cov_matrix_file -m mean_zscores_out_file -e eigenvalues_out_file -q eigenvectors_out_file
7. Simulate zscores from normal distribution with gene covariance matrix.
   - simulate_zscores.py -z input_mean_zscores_file -e input_eigenvalues_file -q input_eigenvectors_file -o out_sim_file -n num_simulations
   - get_sim1000_zscores.py 

