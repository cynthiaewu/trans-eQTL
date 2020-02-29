# trans-eQTL

1. Preprocess the genotype and expression file to get intersecting samples
   - preprocess_genotypefile.py -g input_genotype_file -e input_expression_file -c input_covariates_file -o genotype_output -p expression_output -q covariates_output
2. Filter Genotype file to get only snps in coding regions
   - getCodingRegions.py -i input_preprocess_genotype_file -o output_file -c chr_number
3. Run Matrix eQTL
   - gene-SNP_pairs.R input_coding_genotype_file input_intersect_expression_file pca_file output
4. Get the zscores (T-stat) and pvalues in snp by gene format from Matrix eQTL output
   - get_values.py -i input_matrixeqtl_file -p out_pvalues -z out_zscores
   
   To merge all zscores file into one file:
   
   cp gene_snp_zscores_chr1.csv all_chr_zscores.csv
   
   for CHR in {2..22}; do tail -n +2 -q Nerve-Tibial/chr$CHR/gene_snp_zscores_chr$CHR.csv  >> all_chr_zscores.csv; done
5. Calculate CPMA values for each snp from pvalues file
   - calculate_cpma.py -i input_pvalues_file -o cpma_output
  
    To check the distribution of cpma values follows the chi distribution of df=1, we generated random pvalues from a normal distribution and calculated cpma values for these
   - generate_random_pvalue.py 
6. Calculate the gene covariance matrix and mean zscores for genes with zscores file
   - calculate_cov_meanzscores.py -i input zscores_file -c cov_out -m mzscores_out
7. Calculate the eigendecomposition and mean zscores
   - calculate_edecomposition.py -c input_cov_matrix_file -e eigenvalues_output -q eigenvectors_output
8. Simulate cpma values from normal distribution with gene covariance matrix and mean zscores
   - simulate_zscores.py -z input_mean_zscores_file -e input_eigenvalues_file -q input_eigenvectors_file -o sim_output -n num_simulations
   
   Perform simulations in chunks of 5000
   
   - simulate_cpma_chunks.py -z input_mean_zscores_file -e input_eigenvalues_file -q input_eigenvectors_file -o sim_output -n num_simulations
9. Compare simulated cpma with observed cpma to get an empirical pvalue for each snp
   - calculate_empirical_pvalue.py -s simulated_cpma -o observed_cpma -e empirical_pvalues_output
10. Get qvalues from empirical pvalues
11. Get consequence annotation for every snp.

zcat GTEx_Analysis_20150112_WholeGenomeSeq_VarSitesAnnot.vcf.gz| tail -n +146 | awk -F '\t' '{print $1, $2, $8}' | sed -E 's/(.+) (.+) (.+);CSQ=([^|]+)\|([^|]+)\|(.+)/\1 \2 \5/' > /storage/cynthiawu/trans_eQTL/Gtex_snp-consq.txt
