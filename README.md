# trans-eQTL
## Files
Genotypes:
/storage/mgymrek/gtex-estrs/revision/genotypes/

Expression:
/storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/Review_Rerun

Covariates:
/storage/szfeupe/Runs/GTEx_estr/gtex650.pca

/storage/cynthiawu/trans_eQTL/gtex650_transpose.pca    

## Preprocessing Steps

a. Get gene annotations, keep only protein coding genes

tail -n +7 /storage/resources/dbase/human/gene_annotations/gencode.v19.genes.v7.patched_contigs.gtf| awk -F '\t' '{print $9}' | cut -d ';' -f 1,3,5| cut -d '"' -f 2,4,6 | sed 's/"/\t/g' > /storage/cynthiawu/trans_eQTL/gene_annotations_gencode.v19.genes.v7.csv

b. Get allele count, allele frequency, HWE annotations for each snp

bcftools query -i 'AF[0]>0.01 && AC[0]>2 && HWP>0.001' -f'%CHROM\t%POS\t%AF{1}\t%AC{1}\t%HWP\n' GTEx_Analysis_20150112_WholeGenomeSeq_VarSitesAnnot.vcf.gz -o /storage/cynthiawu/trans_eQTL/GTex_filteredsnps.txt

format_filtersnps.py -i input_bcftool_output -o output_snp_pos_file

```
zcat /storage/resources/datasets/gtex/53844/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v6.p1.c1.GRU/GenotypeFiles/phg000520.v2.GTEx_MidPoint_WGS_SNP_CNV.genotype-calls-vcf.c1/GTEx_Analysis_20150112_WholeGenomeSeq_VarSitesAnnot.vcf.gz| tail -n +146 | awk -F '\t' '{print $1, $2, $8}' | sed -E 's/(.+) (.+) AC=(.+);\AF=(.+);AN=(.+);HWP=(.+);In(.+)/\1 \2 \3 \4 \6/' > /storage/cynthiawu/trans_eQTL/GTEx_snp_AC_AF_HWE.txt

sed -i "1s/.*/chr pos AC AF HWE/" GTEx_snp_AC_AF_HWE.txt
```

c. Get SUBJID, SEX, AGE, TRISCHD, DTHHRDY as covariates

 zcat phs000424.v7.pht002742.v7.p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz| tail -n +11 | awk -F '\t' '{print $2, $4, $5, $15, $35}' OFS='\t' > /storage/cynthiawu/trans_eQTL/gtex_covariates.txt
 
 ## Filtering Scripts
 
 ```
 #combined both steps in filter_genes.py
 Get protein coding genes. 
 - filter_genes.py -i input_file -o out_filtered_file
 
 Get genes with rpkm > 0.1 in at least 10 samples.
 - filter_genes_rpkm.py -i input_file -o out_filtered_file
 
 Get snps with SNPs with MAF >= 0.1%, minor allele count >= 3, HWE >= 0.01
 - filter_snps.py -i input_file -o out_filtered_file
 - getCodingRegions.py -i input_preprocess_genotype_file -o output_file -c chr_number
 
 python filter_snps.py -i ../Nerve-Tibial/chr1/GTExNormalizedSNPGenotypes_chr1_samplename_inter_coding.table -o ../Nerve-Tibial/chr1/GTExNormalizedSNPGenotypes_chr1_samplename_inter_coding_filtered.table
 ```
 
 Get the European sample for the genotype and expression files
 - get_european_samples.py -i input_file -e european_samples_file -t type_file (0 for genotype, 1 for expression) -o data_euro_output
 
 Get the covariates for the European samples
 - get_european_samples_covariates.py -i input_covariate_file -e european_samples_file -o cov_euro_output
 
 ## Steps
1. Intersect expression and genotype files to get intersecting samples. Keep only the intersecting columns of the expression and covariate file. 
   - find_sample_intersect.py -g input_genotype_file -e input_expression_file -c input_covariates_file -i intersect_output -p expression_output -q covariates_output

    Preprocess the genotype files for each chr to get intersecting samples
   - preprocess_genotypefile.py -g input_genotype_file -i intersect_file -o genotype_output 
   
2. Filtering Steps

    Filter Genotype file to get only snps in coding regions. Get snps with SNPs with MAF > 0.01, minor allele count >= 3, HWE > 0.001. snp list = /storage/cynthiawu/trans_eQTL/GTex_filteredsnps_pos.txt
   - filter_snps.py -i input_preprocess_genotype_file -s filtered_snp_list -o output_file -c chr_number
   
    Filter Expression file. Get protein coding genes and genes with rpkm > 0.1 in at least 10 samples.
   - filter_genes.py -i input_file -o out_filtered_file
 
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
