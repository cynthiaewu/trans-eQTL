import sys
import subprocess
import os
import argparse

def cpma_pipeline(genotype, expression, covariates, scripts_folder, output):
    #genotype = f'{input_folder}/genotype.csv'
    #expression = f'{input_folder}/expression.csv'
    #genotype = f'{input_folder}/genotype_permuteSNP0.csv'
    #expression = f'{input_folder}/expression.csv'
    cpma_folder = os.path.join(output, 'CPMA')
    cpmax_folder = os.path.join(output, 'CPMAx')
    if not os.path.isdir(cpma_folder):
        os.mkdir(cpma_folder)
    if not os.path.isdir(cpmax_folder):
        os.mkdir(cpmax_folder)
    eqtl_file = f'{cpma_folder}/gene-snp-eqtl'

    matrix_cmd = f'Rscript /storage/polo/GTEx_v8_matrix_eQTL/Scripts/Run_Matrix_eQTL_PC_PEER_quantile_norm.r {genotype} {expression} {covariates} {eqtl_file}'.split(' ')

    #subprocess.call(matrix_cmd)
    print(f'Finished matrix eQTL')
    # Obtain zscores and pvalues in a snp by gene matrix format from matrix eQTL output
    pvalue_file = f'{eqtl_file}_pvalue'
    zscore_file = f'{eqtl_file}_zscore'

    values_cmd = f'python {scripts_folder}/CPMA/get_values.py -i {eqtl_file} -z {zscore_file}'.split(' ')
    values = subprocess.Popen(values_cmd).wait()
    print('Finished getting pvalues and zscores')
    
    
    pvalue_file = f'{eqtl_file}_pvalue_converted'
    zscore_file = f'{eqtl_file}_zscore'

    values_cmd = f'python {scripts_folder}/CPMA/tstat_to_pvalue.py -i {zscore_file} -p {pvalue_file}'.split(' ')
    values = subprocess.Popen(values_cmd).wait()
    print('Finished converting zscores to pvalues')
    
    #Calculate cpma values for all chr with matrix eqtl pvalues output
    cpma_file = f'{eqtl_file}_cpma_converted'
    
    cpma_cmd = f'python {scripts_folder}/CPMA/calculate_cpma.py -i {pvalue_file} -o {cpma_file}'.split(' ')
    cpma = subprocess.Popen(cpma_cmd).wait()
    print('Finished calculating cpma')
    ''' 
    #Calculate cpma values for all chr with matrix eqtl pvalues output
    cpma_file = f'{eqtl_file}_cpma-mix_converted'
    
    #pvalue_file = f'{eqtl_file}_pvalue_converted'
    #for i in range(1, 23):
    cpma_cmd = f'python ../mixtureModel/calculate_mixturemodel.py -i {output}/chr{i}/{pvalue_file}_chr{i} -o {output}/chr{i}/{cpma_file}_chr{i}'.split(' ')
    cpma = subprocess.Popen(cpma_cmd).wait()
    print(f'Finished calculating cpma-mix for chr{i}')
    print('Finished calculating cpma-mix for all chr')
    
    cpma_file = f'{eqtl_file}_cpma_converted'
   
    num_sim = 500000
    #Compare simulated cpma with observed cpma to get an empirical pvalue for each snp
    compare_cmd = []
    #for i in range(1, 23):
    for i in [1]:
        sim_file = f'{output}/chr{i}/{eqtl_file}_sim{num_sim}_cpma_chr{i}_converted'
        empirical_file = f'{output}/chr{i}/{eqtl_file}_empiricalpvalues_chr{i}_converted'
        compare_cmd.append(f'python calculate_empirical_pvalue.py -s {sim_file} -o {output}/chr{i}/{cpma_file}_chr{i} -e {empirical_file}'.split(' '))
        if len(compare_cmd) > 0 or (i==22):
            compare_procs = [ subprocess.Popen(i) for i in compare_cmd]
            print(compare_procs)
            for p in compare_procs:
                p.wait()
            compare_cmd = []
    print('Finished calculating empirical pvalues from cpma for all chr')
    
    #Merge all empirical pvalues of all chr to one file
    all_epvalues = f'{output}/{eqtl_file}_empiricalpvalues_allchr_converted'
    cp_cmd = f'cp {output}/chr1/{eqtl_file}_empiricalpvalues_chr1_converted {all_epvalues}'.split(' ')
    cp = subprocess.Popen(cp_cmd).wait()
    for i in range(2, 23):
        append_cmd = f'tail -n +2 -q {output}/chr{i}/{eqtl_file}_empiricalpvalues_chr{i}_converted  >> {all_epvalues}'
        subprocess.call(append_cmd, shell=True)
    print('Finished merging all chr empirical pvalues to one file')
    
    for chrm, error in chr_error.items():
        if error != 0:
            print(f'chr{chrm}: Did not finish, Errorcode: {error}')
    '''

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-g', '--genotype', type=str, help='Input genotype file')
    parser.add_argument('-e', '--expression', type=str, help='Input expression file')
    parser.add_argument('-c', '--covariates', type=str, help='Input covariates file')
    parser.add_argument('-p', '--scripts_folder', type=str, help='Input folder with scripts')
    parser.add_argument('-o', '--output', type=str, help='Destination folder for output')
    params = parser.parse_args()

    cpma_pipeline(params.genotype, params.expression, params.covariates, params.scripts_folder, params.output)


if __name__ == '__main__':
    main()
