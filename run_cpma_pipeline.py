import sys
import subprocess
import os
import argparse

def cpma_pipeline(genotype, expression, tissue, covariate, european, output):
    genotype_file1 = genotype + '/GTExNormalizedSNPGenotypes_chr1.table.gz'
    intersect_file = f'{output}/intersect_samples.txt'
    expression_inter = f'{output}/Clean_expression_{tissue}_inter'
    covariates_inter = f'{output}/gtex650_{tissue}_intersect.pca'
    
    #Get the intersecting samples in the expression, genotype, and covariates file
    intersect_cmd = f'python find_sample_intersect.py -g {genotype_file1} -e {expression} -c {covariate} -i {intersect_file} -p {expression_inter} -q {covariates_inter}'.split(' ')
    intersect = subprocess.Popen(intersect_cmd).wait()
    print('Intersecting Files Done')

    for i in range(1, 23):
        output_path = os.path.join(output, f'chr{i}')
        if not os.path.isdir(output_path):
            os.mkdir(output_path)

    #Get obly the intersecting samples in the genotype files for all chr
    genotype_inter = f'GTExNormalizedSNPGenotypes_inter'
    for i in range(1, 23):
        preprocess_cmd = f'python preprocess_genotypefile.py -g {genotype}/GTExNormalizedSNPGenotypes_chr{i}.table.gz -i {intersect_file} -o {output}/chr{i}/{genotype_inter}_chr{i}'.split(' ')
        #print(preprocess_cmd)
        preprocess = subprocess.Popen(preprocess_cmd).wait()
        print(f'Finished preprocessing chr{i}')
    print('Preprocessing Genotype Files Done')

    #Filter snps in genotype files for all chr
    genotype_filter = f'{genotype_inter}_filter'
    filtered_snp_list = '/storage/cynthiawu/trans_eQTL/GTex_filteredsnps_pos.txt'
    for i in range(1, 23):
        filtersnp_cmd = f'python filter_snps.py -i {output}/chr{i}/{genotype_inter}_chr{i} -s {filtered_snp_list} -o {output}/chr{i}/{genotype_filter}_chr{i} -c {i}'.split(' ')
        filtersnp = subprocess.Popen(filtersnp_cmd).wait()
        print(f'Finished Filtering chr{i}')
    print('Filtering Snps Done')
    
    #Filter genes in expression file
    expression_filter = f'{expression_inter}_filter'
    filtergene_cmd = f'python filter_genes.py -i {expression_inter} -o {expression_filter}'.split(' ')
    filtergene = subprocess.Popen(filtergene_cmd).wait()
    print('Filtering Genes Done')
  
    #Perform matrix eQTL on all chr to get gene-snp pairs
    eqtl_file = 'gene-snp_eqtls'
    for i in range(1, 23):
        matrix_cmd = f'Rscript gene-SNP_pairs.R {output}/chr{i}/{genotype_filter}_chr{i} {expression_filter} {covariates_inter} {output}/chr{i}/{eqtl_file}_chr{i}'.split(' ')
        matrix = subprocess.Popen(matrix_cmd).wait()
        print(f'Finished matrix eQTL for chr{i}')
    print('Finished matrix eQTL for all chr')

    # Obtain zscores and pvalues in a snp by gene matrix format from matrix eQTL output
    wc_cmd = f'wc -l {expression_filter}'.split(' ')
    wc = subprocess.Popen(wc_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = wc.communicate()
    num_genes = int(stdout.split()[0])-1
    pvalue_file = f'{eqtl_file}_pvalue'
    zscore_file = f'{eqtl_file}_zscore'
    for i in range(1, 23):
        values_cmd = f'python get_values.py -i {output}/chr{i}/{eqtl_file}_chr{i} -n {num_genes} -p {output}/chr{i}/{pvalue_file}_chr{i} -z {output}/chr{i}/{zscore_file}_chr{i}'.split(' ')
        values = subprocess.Popen(values_cmd).wait()
        print(f'Finished getting pvalues and zscores for chr{i}')
    print('Finished getting pvalues and zscores for all chr')

    #Calculate cpma values for all chr with matrix eqtl pvalues output
    cpma_file = f'{eqtl_file}_cpma'
    for i in range(1, 23):
        cpma_cmd = f'python calculate_cpma.py -i {output}/chr{i}/{pvalue_file}_chr{i} -o {output}/chr{i}/{cpma_file}_chr{i}'.split(' ')
        cpma = subprocess.Popen(cpma_cmd).wait()
        print(f'Finished calculating cpma for chr{i}')
    print('Finished calculating cpma for all chr')

    #Combine zscores of all chr in one file
    all_out = f'{output}/{eqtl_file}_allchr'
    all_zscores = f'{all_out}_zscores'
    cp_cmd = f'cp {output}/chr1/{zscore_file}_chr1 {all_zscores}'.split(' ')
    cp = subprocess.Popen(cp_cmd).wait()
    for i in range(2, 23):
        append_cmd = f'tail -n +2 -q {output}/chr{i}/{zscore_file}_chr{i}  >> {all_zscores}'
        subprocess.call(append_cmd, shell=True)
    print('Finished merging all chr zscores to one file')

    #Calculate the gene covariance matrix and mean zscores for genes
    cov_matrix = f'{all_out}_cov'
    mzscores = f'{all_out}_meanzscores'
    cov_cmd = f'python calculate_cov_meanzscores.py -i {all_zscores} -c {cov_matrix} -m {mzscores}'.split(' ')
    cov = subprocess.Popen(cov_cmd).wait()
    print('Finished calculating cov matrix and mean zscores')

    #Calculate the eigendecomposition
    evalues_file = f'{all_out}_evalues'
    evectors_file = f'{all_out}_Q'
    eigen_cmd = f'python calculate_edecomposition.py -c {cov_matrix} -e {evalues_file} -q {evectors_file}'.split(' ')
    eigen = subprocess.Popen(eigen_cmd).wait()
    print('Finished calculating eigendecomposition')

    #Simulate cpma values from normal distribution with gene covariance matrix and mean zscores
    num_sim = 5000000
    sim_file = f'{all_out}_sim{num_sim}_cpma'
    simulate_cmd = f'python simulate_cpma_chunks.py -z {mzscores} -e {evalues_file} -q {evectors_file} -o {sim_file} -n {num_sim}'.split(' ')
    simulate = subprocess.Popen(simulate_cmd).wait()
    print('Finished simulation of cpma values')

    #Compare simulated cpma with observed cpma to get an empirical pvalue for each snp
    empirical_file = f'{eqtl_file}_empiricalpvalues'
    for i in range(1, 23):
        compare_cmd = f'python calculate_empirical_pvalue.py -s {sim_file} -o {output}/chr{i}/{cpma_file}_chr{i} -e {output}/chr{i}/{empirical_file}_chr{i}'.split(' ')
        compare = subprocess.Popen(compare_cmd).wait()
    print('Finished calculating empirical pvalues from cpma for all chr')

    #Merge all empirical pvalues of all chr to one file
    all_epvalues = f'{all_out}_empiricalpvalues'
    cp_cmd = f'cp {output}/chr1/{empirical_file}_chr1 {all_epvalues}'.split(' ')
    cp = subprocess.Popen(cp_cmd).wait()
    for i in range(2, 23):
        append_cmd = f'tail -n +2 -q {output}/chr{i}/{empirical_file}_chr{i}  >> {all_epvalues}'
        subprocess.call(append_cmd, shell=True)
    print('Finished merging all chr empirical pvalues to one file')
    '''
    finalgenofile = genotype_filter
    finalexpfile = expression_filter
    if european == 1:
        eurosamp = '/storage/cynthiawu/trans_eQTL/European_samples.txt'
        finalgenofile = f'{genotype_filter}_european'
        for i in range(1, 3):
            eurogeno_cmd = f'python get_european_samples.py -i {output}/chr{i}/{genotype_filter}_chr{i} -e {eurosamp} -t 0 -o {output}/chr{i}/{finalgenofile}_chr{i}'.split(' ')
            print(eurogeno_cmd)
            eurogeno = subprocess.Popen(eurogeno_cmd).wait()
            print(f'Obtained European samples for chr{i}')
        finalexpfile = f'{expression_filter}_european'
        euroexp_cmd = f'python get_european_samples.py -i {expression_filter} -e {eurosamp} -t 1 -o {finalexpfile}'.split(' ')
        euroexp = subprocess.Popen(euroexp_cmd).wait()
        print('Obtained European samples for expression data')
        eurocovar_cmd = f'python get_european_samples_covariates.py -i {covariates_inter} -e {eurosamp} -o {output}/gtex650_eurosamp_covariates.pca'.split(' ')
        eurocovar = subprocess.Popen(eurocovar_cmd).wait()
        print('Obtained European samples for covariates')
    '''
    '''
    readcounts = os.path.join(output_path, sample + "_readcounts.txt")
    subprocess.call(["python", "prepare_readcounts.py", "-i", readcounts_path, "-o", readcounts])
    '''

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-g', '--genotype', type=str, help='Input folder with genotype files')
    parser.add_argument('-e', '--expression', type=str, help='Input expression file')
    parser.add_argument('-t', '--tissue', type=str, help='Tissue type')
    parser.add_argument('-c', '--covariate', type=str, help='Input covariate file')
    parser.add_argument('-s', '--european', type=int, help='0: Include all samples, 1: Include only European samples')
    parser.add_argument('-o', '--output', type=str, help='Destination folder for output')
    input_values = parser.parse_args()

    cpma_pipeline(input_values.genotype, input_values.expression, input_values.tissue, input_values.covariate, input_values.european, input_values.output)


if __name__ == '__main__':
    main()
