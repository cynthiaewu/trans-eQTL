# trans-eQTL

1. Run Matrix eQTL
  - gene-SNP_pairs.R 
2. Get the zscores (T-stat) and pvalues from Matrix eQTL output
  - get_pvalues.py
3. Calculate CPMA values for each snp from pvalues file
  - get_cpma.py
  
    To check the distribution of cpma values follows the chi distribution of df=1, we generated random pvalues from a normal distribution and calculated cpma values for these
  - generate_random_pvalue.py 
4. Calculate the gene covariance matrix with zscores file
  - get_cov_matrix.py 
5. Simulate zscores from normal distribution with gene covariance matrix.
  - get_sim_zscores.py
  - get_sim1000_zscores.py 

