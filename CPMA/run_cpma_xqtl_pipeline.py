import sys
import subprocess
import os
import argparse

def cpmaxqtl_pipeline(input_folder, scripts_folder, topx, method, matrixeqtl, genecorr):
    genotype = f'{input_folder}/genotype.csv'
    expression = f'{input_folder}/expression.csv'
    matrixeqtl_folder = os.path.join(input_folder, 'Matrixeqtl')
    cpma_folder = os.path.join(input_folder, 'CPMA')
    xqtl_folder = os.path.join(input_folder, 'x-QTL')
    if not os.path.isdir(matrixeqtl_folder):
        os.mkdir(matrixeqtl_folder)
    if not os.path.isdir(cpma_folder):
        os.mkdir(cpma_folder)
    if not os.path.isdir(xqtl_folder):
        os.mkdir(xqtl_folder)
    eqtl_file = f'{matrixeqtl_folder}/gene-snp-eqtl'    
    #eqtl_file = f'{cpma_folder}/yeast_matrixeqtl_PEER'    
    
    if matrixeqtl: 
        #Perform matrix eQTL to get gene-snp pairs
        matrix_cmd = f'Rscript {scripts_folder}/MatrixeQTL/Run_Matrix_eQTL_PC_PEER_quantile_norm.r -g {genotype} -e {expression} -o {eqtl_file}'.split(' ')
        subprocess.call(matrix_cmd)
        print(f'Finished matrix eQTL, {input_folder}')   
    
    #Obtain zscores and pvalues in a snp by gene matrix format from matrix eQTL output    
    pvalue_file = f'{eqtl_file}_pvalue_converted.gz'
    zscore_file = f'{eqtl_file}_zscore.gz'
    #pvalue_file = f'{input_folder}/yeast_pvals.csv'
    #zscore_file = f'{input_folder}/yeast_zscores.csv'
     
    #values_cmd = f'python /storage/cynthiawu/trans_eQTL/Scripts/CPMA/get_values.py -i {eqtl_file} -n {num_genes} -p {pvalue_file} -z {zscore_file}'.split(' ')
    values_cmd = f'python {scripts_folder}/CPMA/get_values.py -i {eqtl_file} -z {zscore_file}'.split(' ')
    values = subprocess.Popen(values_cmd).wait()
    print(f'Finished getting pvalues and zscores, {input_folder}')
    
    convert_cmd = f'python {scripts_folder}/CPMA/tstat_to_pvalue.py -i {zscore_file} -p {pvalue_file}'.split(' ')
    convert = subprocess.Popen(convert_cmd).wait()
    
    if method==0 or method==1:
        #Calculate cpma values with matrix eqtl pvalues output
        cpma_file = f'{cpma_folder}/gene-snp-eqtl_cpma_topx_{topx}_converted'
        
        #cpma_cmd = f'python /storage/cynthiawu/trans_eQTL/Scripts/CPMA/calculate_cpma_topx.py -i {pvalue_file} -x {topx} -o {cpma_file}'.split(' ')
        cpma_cmd = f'python {scripts_folder}/CPMA/calculate_cpma_topx.py -i {pvalue_file} -x {topx} -o {cpma_file}'.split(' ')
        cpma = subprocess.Popen(cpma_cmd).wait()
        print(f'Finished calculating cpma, {input_folder}')
    
    if method==0 or method==2:
        xqtl_file = f'{xqtl_folder}/gene-snp-eqtl_xqtl'
        xqtl_cmd = f'python {scripts_folder}/x-QTL/calculate_mixturemodel.py -i {pvalue_file} -o {xqtl_file}'.split(' ')
        xqtl = subprocess.Popen(xqtl_cmd).wait()
        print(f'Finished calculating x-QTL')
    
    if genecorr:
        mzscores = f'{eqtl_file}_meanzscores.gz'
        evalues_file = f'{eqtl_file}_evalues.gz'
        evectors_file = f'{eqtl_file}_Q.gz'
        
        cov_cmd = (f'python {scripts_folder}/CPMA/calculate_cov_meanzscores_edecomposition.py -i {zscore_file} -m {mzscores} -e {evalues_file} -q {evectors_file}'.split(' '))
        cov = subprocess.Popen(cov_cmd).wait()
        print('Finished calculating mean zscores and eigendecomposition')
        
        #Simulate cpma values from normal distribution with gene covariance matrix and mean zscores
        num_sim = 500000
        sim_file = f'{eqtl_file}_sim{num_sim}_cpma.gz'
        sim_mix_file = f'{eqtl_file}_sim{num_sim}_xqtl.gz'
        simulate_cmd = f'python {scripts_folder}/x-QTL/simulate_cpma-mix_chunks.py -z {mzscores} -e {evalues_file} -q {evectors_file} -o {sim_file} -m {sim_mix_file} -n {num_sim}'.split(' ')
        #simulate_cmd = f'python {scripts_folder}/CPMA/simulate_cpma_chunks.py -z {mzscores} -e {evalues_file} -q {evectors_file} -o {sim_file} -n {num_sim}'.split(' ')
        simulate = subprocess.Popen(simulate_cmd).wait()
        print('Finished simulation of cpma values')
        
        if method==0 or method==1:
            cpma_file = f'{cpma_folder}/gene-snp-eqtl_cpma_topx_{topx}_converted'
            num_sim = 500000
            sim_file = f'{eqtl_file}_sim{num_sim}_cpma.gz'
            #Compare simulated cpma with observed cpma to get an empirical pvalue for each snp
            empirical_file = f'{eqtl_file}_empiricalpvalues_topx_{topx}_converted'
            compare_cmd = f'python {scripts_folder}/CPMA/calculate_empirical_pvalue.py -s {sim_file} -o {cpma_file} -e {empirical_file}'.split(' ')
            compare = subprocess.Popen(compare_cmd).wait()
            print('Finished calculating empirical pvalues for cpma from simulated empirical null distribution, gene correlation')
    else:
        if method==0 or method==1: 
            cpma_file = f'{cpma_folder}/gene-snp-eqtl_cpma_topx_{topx}_converted'
            empirical_file = f'{eqtl_file}_chidist_topx_{topx}_converted'
            emp_cmd = (f'python {scripts_folder}/Simulator/compute_pvalue_chidist_nocorr.py -i {cpma_file} -o {empirical_file}'.split(' '))
            emp = subprocess.Popen(emp_cmd).wait()
            print('Finished calculating empirical pvalues for cpma from chi distribution, no gene correlation')


def iterate_folders(folder, scripts_folder, topx, method, iterations, matrixeqtl, genecorr):
    for i in range(iterations):
        input_folder = f'{folder}/Simulation_{i}'
        print(f'Starting eqtl for Simulation {i}, {folder}')
        cpmaxqtl_pipeline(input_folder, scripts_folder, topx, method, matrixeqtl, genecorr)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--folder', type=str, help='Input folder with simulated expression and genotype files')
    parser.add_argument('-p', '--scripts_folder', type=str, help='Input folder with scripts')
    parser.add_argument('-x', '--topx', default=1.0, type=float, help='Top x of genes to be used in cpma calculation (0 to 1)')
    parser.add_argument('-m', '--method', default=0, type=int, help='CPMA=1, x-QTL=2, both=0')
    parser.add_argument('-i', '--iterations', type=int, help='# of simulated folders to run cpma pipeline on')
    parser.add_argument('-q', '--matrixeqtl', action='store_true', help='Run Matrix eQTL')
    parser.add_argument('-g', '--genecorr', action='store_true', help='Run eigendecomposition to account for gene correlation')
    params = parser.parse_args()
    
    iterate_folders(folder=params.folder, scripts_folder=params.scripts_folder, topx=params.topx, method=params.method, iterations=params.iterations, matrixeqtl=params.matrixeqtl, genecorr=params.genecorr)


if __name__ == '__main__':
    main()
