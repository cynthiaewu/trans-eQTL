import sys
import subprocess
import os
import argparse

def cpmax_pipeline(input_folder, scripts_folder, topx):
    genotype = f'{input_folder}/genotype.csv'
    expression = f'{input_folder}/expression.csv'
    cpma_folder = os.path.join(input_folder, 'CPMA')
    cpmax_folder = os.path.join(input_folder, 'CPMAx')
    if not os.path.isdir(cpma_folder):
        os.mkdir(cpma_folder)
    if not os.path.isdir(cpmax_folder):
        os.mkdir(cpmax_folder)
    eqtl_file = f'{cpma_folder}/gene-snp-eqtl'
    
    
    #Perform matrix eQTL to get gene-snp pairs
    #matrix_cmd = f'Rscript /storage/cynthiawu/trans_eQTL/Scripts/MatrixeQTL/gene-SNP_pairs.R -g {genotype} -e {expression} -o {eqtl_file}'.split(' ')
    matrix_cmd = f'Rscript {scripts_folder}/MatrixeQTL/gene-SNP_pairs.R -g {genotype} -e {expression} -o {eqtl_file}'.split(' ')
    subprocess.call(matrix_cmd)
    print(f'Finished matrix eQTL, {input_folder}')
    
    
    #Obtain zscores and pvalues in a snp by gene matrix format from matrix eQTL output    
    pvalue_file = f'{eqtl_file}_pvalue.gz'
    zscore_file = f'{eqtl_file}_zscore.gz'
    
    '''
    #values_cmd = f'python /storage/cynthiawu/trans_eQTL/Scripts/CPMA/get_values.py -i {eqtl_file} -n {num_genes} -p {pvalue_file} -z {zscore_file}'.split(' ')
    values_cmd = f'python {scripts_folder}/CPMA/get_values.py -i {eqtl_file} -p {pvalue_file} -z {zscore_file}'.split(' ')
    values = subprocess.Popen(values_cmd).wait()
    print(f'Finished getting pvalues and zscores, {input_folder}')
    '''

    
    #Calculate cpma values with matrix eqtl pvalues output
    cpma_file = f'{cpmax_folder}/gene-snp-eqtl_cpma_topx_{topx}'
    '''
    #cpma_cmd = f'python /storage/cynthiawu/trans_eQTL/Scripts/CPMA/calculate_cpma_topx.py -i {pvalue_file} -x {topx} -o {cpma_file}'.split(' ')
    cpma_cmd = f'python {scripts_folder}/CPMA/calculate_cpma_topx.py -i {pvalue_file} -x {topx} -o {cpma_file}'.split(' ')
    cpma = subprocess.Popen(cpma_cmd).wait()
    print(f'Finished calculating cpma, {input_folder}')
    '''
   
    #Calculate the gene covariance matrix and mean zscores for genes
    cov_matrix = f'{eqtl_file}_cov.gz'
    
    mzscores = f'{eqtl_file}_meanzscores.gz'
    evalues_file = f'{eqtl_file}_evalues.gz'
    evectors_file = f'{eqtl_file}_Q.gz'
    
    cov_cmd = f'python {scripts_folder}/CPMA/calculate_cov_meanzscores_edecomposition.py -i {zscore_file} -m {mzscores} -e {evalues_file} -q {evectors_file}'.split(' ')
    cov = subprocess.Popen(cov_cmd).wait()
    print('Finished calculating mean zscores and eigendecomposition')

    ''' 
    #Calculate the eigendecomposition
    evalues_file = f'{eqtl_file}_evalues.gz'
    evectors_file = f'{eqtl_file}_Q.gz'
    
    eigen_cmd = f'python {scripts_folder}/CPMA/calculate_edecomposition.py -c {cov_matrix} -e {evalues_file} -q {evectors_file}'.split(' ')
    eigen = subprocess.Popen(eigen_cmd).wait()
    print('Finished calculating eigendecomposition')
    
    rm_cmd = f'rm {cov_matrix}'.split(' ')
    rm = subprocess.Popen(rm_cmd).wait()
    '''
    '''
    #Simulate cpma values from normal distribution with gene covariance matrix and mean zscores
    num_sim = 500000
    sim_file = f'{eqtl_file}_sim{num_sim}_cpma.gz'
    simulate_cmd = f'python -u {scripts_folder}/CPMA/simulate_cpma_chunks.py -z {mzscores} -e {evalues_file} -q {evectors_file} -o {sim_file} -n {num_sim}'.split(' ')
    simulate = subprocess.Popen(simulate_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    for line in simulate.stdout:
        sys.stdout.write(line)
    simulate.wait()
    print('Finished simulation of cpma values')

    #Compare simulated cpma with observed cpma to get an empirical pvalue for each snp
    empirical_file = f'{eqtl_file}_empiricalpvalues_topx_{topx}'
    compare_cmd = f'python {scripts_folder}/CPMA/calculate_empirical_pvalue.py -s {sim_file} -o {cpma_file} -e {empirical_file}'.split(' ')
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
