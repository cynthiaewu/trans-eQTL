import sys
import subprocess
import os
import argparse

def cpmax_pipeline(input_folder, scripts_folder, topx):
    genotype = f'{input_folder}/genotype.csv'
    expression = f'{input_folder}/expression_PCs.csv'
    #covariates = f'{input_folder}/covariates_geno.csv'
    #cpma_folder = os.path.join(input_folder, 'CPMA')
    #cpmax_folder = os.path.join(input_folder, 'CPMAx')
    PCs_folder = os.path.join(input_folder, 'expressionPCs')
    if not os.path.isdir(PCs_folder):
        os.mkdir(PCs_folder)
    #eqtl_file = f'{PCs_folder}/gene-snp-eqtl_PCs'
    eqtl_file = f'{PCs_folder}/yeast_matrixeqtl_PCs_PEER'
    '''
    #Perform matrix eQTL to get gene-snp pairs
    #matrix_cmd = f'Rscript /storage/cynthiawu/trans_eQTL/Scripts/MatrixeQTL/gene-SNP_pairs.R -g {genotype} -e {expression} -o {eqtl_file}'.split(' ')
    matrix_cmd = f'Rscript {scripts_folder}/MatrixeQTL/gene-SNP_pairs.R -g {genotype} -e {expression} -o {eqtl_file}'.split(' ')
    #matrix_cmd = f'Rscript {scripts_folder}/MatrixeQTL/gene-SNP_pairs.R -g {genotype} -e {expression} -c {covariates} -o {eqtl_file}'.split(' ')
    subprocess.call(matrix_cmd)
    print(f'Finished matrix eQTL, {input_folder}')
    '''
    
    #Obtain zscores and pvalues in a snp by gene matrix format from matrix eQTL output
    pvalue_file = f'{eqtl_file}_pvalue'
    zscore_file = f'{eqtl_file}_zscore'
   
    #values_cmd = f'python /storage/cynthiawu/trans_eQTL/Scripts/CPMA/get_values.py -i {eqtl_file} -n {num_genes} -p {pvalue_file} -z {zscore_file}'.split(' ')
    values_cmd = f'python {scripts_folder}/CPMA/get_values.py -i {eqtl_file} -z {zscore_file}'.split(' ')
    values = subprocess.Popen(values_cmd).wait()
    values_cmd = f'python {scripts_folder}/CPMA/tstat_to_pvalue.py -i {zscore_file} -p {pvalue_file}'.split(' ')   
    values = subprocess.Popen(values_cmd).wait()
    print(f'Finished getting pvalues and zscores, {input_folder}')

    #Calculate cpma values with matrix eqtl pvalues output
    cpma_file = f'{eqtl_file}_cpma_topx_{topx}'
    #cpma_cmd = f'python /storage/cynthiawu/trans_eQTL/Scripts/CPMA/calculate_cpma_topx.py -i {pvalue_file} -x {topx} -o {cpma_file}'.split(' ')
    cpma_cmd = f'python {scripts_folder}/CPMA/calculate_cpma_topx.py -i {pvalue_file} -x {topx} -o {cpma_file}'.split(' ')
    cpma = subprocess.Popen(cpma_cmd).wait()
    print(f'Finished calculating cpma, {input_folder}')
    
    '''
    #Calculate the gene covariance matrix and mean zscores for genes
    cov_matrix = f'{eqtl_file}_cov'
    mzscores = f'{eqtl_file}_meanzscores'
    cov_cmd = f'python calculate_cov_meanzscores.py -i {zscore_file} -c {cov_matrix} -m {mzscores}'.split(' ')
    cov = subprocess.Popen(cov_cmd).wait()
    print('Finished calculating cov matrix and mean zscores')

    #Calculate the eigendecomposition
    evalues_file = f'{eqtl_file}_evalues'
    evectors_file = f'{eqtl_file}_Q'
    eigen_cmd = f'python calculate_edecomposition.py -c {cov_matrix} -e {evalues_file} -q {evectors_file}'.split(' ')
    eigen = subprocess.Popen(eigen_cmd).wait()
    print('Finished calculating eigendecomposition')

    #Simulate cpma values from normal distribution with gene covariance matrix and mean zscores
    num_sim = 1000000
    sim_file = f'{eqtl_file}_sim{num_sim}_cpma'
    simulate_cmd = f'python simulate_cpma_chunks.py -z {mzscores} -e {evalues_file} -q {evectors_file} -o {sim_file} -n {num_sim}'.split(' ')
    simulate = subprocess.Popen(simulate_cmd).wait()
    print('Finished simulation of cpma values')

    #Compare simulated cpma with observed cpma to get an empirical pvalue for each snp
    empirical_file = f'{eqtl_file}_empiricalpvalues'
    compare_cmd = f'python calculate_empirical_pvalue.py -s {sim_file} -o {cpma_file} -e {empirical_file}'.split(' ')
    compare = subprocess.Popen(compare_cmd).wait()
    print('Finished calculating empirical pvalues from cpma')
    '''


def iterate_folders(folder, scripts_folder, topx, iterations):
    for i in range(iterations):
        input_folder = f'{folder}/Simulation_{i}'
        print(f'Starting CPMA for Simulation {i}, {folder}')
        cpmax_pipeline(input_folder, scripts_folder, topx)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--folder', type=str, help='Input folder with simulated expression and genotype files')
    parser.add_argument('-p', '--scripts_folder', type=str, help='Input folder with scripts')
    parser.add_argument('-x', '--topx', default=0.1, type=float, help='Top x of genes to be used in cpma calculation (0 to 1)')
    parser.add_argument('-i', '--iterations', type=int, help='# of simulated folders to run cpma pipeline on')
    params = parser.parse_args()

    iterate_folders(folder=params.folder, scripts_folder=params.scripts_folder, topx=params.topx, iterations=params.iterations)


if __name__ == '__main__':
    main()
