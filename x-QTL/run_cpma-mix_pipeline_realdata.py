import sys
import subprocess
import os
import argparse

def cpma_pipeline(input_folder, output):
    eqtl_file = f'gene-snp_eqtls'
    
    if not os.path.isdir(output):
        os.mkdir(output)
    
    for i in range(1, 23):
    #for i in [1]:
        chr_folder = os.path.join(output, f'chr{i}')
        if not os.path.isdir(chr_folder):
            os.mkdir(chr_folder)
    
   
    # Obtain zscores and pvalues in a snp by gene matrix format from matrix eQTL output
    for i in range(1, 23):    
    #for i in [1]:
        pvalue_file = f'{eqtl_file}_pvalue'
        zscore_file = f'{eqtl_file}_zscore'
    
        #values_cmd = f'python get_values.py -i {output}/chr{i}/matrixeqtl/matrixeQTL_results_chr{i}.gz -z {output}/chr{i}/CPMA/{zscore_file}_chr{i}'.split(' ')
        values_cmd = f'python get_values.py -i {output}/chr{i}/matrixeqtl_results_ciseqtls_chr{i}.gz -z {output}/chr{i}/CPMA/{zscore_file}_chr{i}'.split(' ')
        values = subprocess.Popen(values_cmd).wait()
        print(f'Finished getting pvalues and zscores for chr{i}')
    print('Finished getting pvalues and zscores for all chr')
    
    
    for i in range(1, 23):
        pvalue_file = f'{eqtl_file}_pvalue_converted'
        zscore_file = f'{eqtl_file}_zscore'
    
        values_cmd = f'python tstat_to_pvalue.py -i {output}/chr{i}/CPMA/{zscore_file}_chr{i} -p {output}/chr{i}/CPMA/{pvalue_file}_chr{i}'.split(' ')
        values = subprocess.Popen(values_cmd).wait()
        print(f'Finished converting zscores to pvalues for chr{i}')
    print('Finished converting zscores to pvalues for all chr')
    
    #Calculate cpma values for all chr with matrix eqtl pvalues output
    cpma_file = f'{eqtl_file}_cpma_converted'
    
    for i in range(1, 23):
        cpma_cmd = f'python calculate_cpma.py -i {output}/chr{i}/CPMA/{pvalue_file}_chr{i} -o {output}/chr{i}/CPMA/{cpma_file}_chr{i}'.split(' ')
        cpma = subprocess.Popen(cpma_cmd).wait()
        print(f'Finished calculating cpma for chr{i}')
    print('Finished calculating cpma for all chr')
    
    #Calculate cpma values for all chr with matrix eqtl pvalues output
    pvalue_file = f'{eqtl_file}_pvalue_converted'
    cpma_file = f'{eqtl_file}_cpma-mix_converted'
    
    #pvalue_file = f'{eqtl_file}_pvalue_converted'
    for i in range(1, 23):
        cpma_cmd = f'python ../mixtureModel/calculate_mixturemodel.py -i {output}/chr{i}/CPMA/{pvalue_file}_chr{i} -o {output}/chr{i}/CPMA/{cpma_file}_chr{i}'.split(' ')
        cpma = subprocess.Popen(cpma_cmd).wait()
        print(f'Finished calculating cpma-mix for chr{i}')
    print('Finished calculating cpma-mix for all chr')
    
    '''
    cov_cmd = []
    for i in range(1, 23):
        print('starting eigendecomp')
        #cov_matrix = f'{output}/chr{i}/{eqtl_file}_cov_chr{i}'
        mzscores = f'{output}/chr{i}/CPMA/{eqtl_file}_meanzscores_chr{i}.gz'
        evalues_file = f'{output}/chr{i}/CPMA/{eqtl_file}_evalues_chr{i}.gz'
        evectors_file = f'{output}/chr{i}/CPMA/{eqtl_file}_Q_chr{i}.gz'
        cov_cmd.append(f'python calculate_cov_meanzscores_edecomposition.py -i {output}/chr{i}/CPMA/{eqtl_file}_zscore_chr{i} -m {mzscores} -e {evalues_file} -q {evectors_file}'.split(' '))

        if len(cov_cmd) > 1 or (i==22):
            cov_procs = [ subprocess.Popen(i) for i in cov_cmd]
            print(cov_procs)
            for p in cov_procs:
                p.wait()
            cov_cmd = []
        #cov = subprocess.Popen(cov_cmd).wait()
    print('Finished calculating mean zscores and eigendecomposition')
    '''
    '''
    #Simulate cpma values from normal distribution with gene covariance matrix and mean zscores
    simulate_cmd = []
    chr_error = {}
    #for i in range(1, 23):
        mzscores = f'{output}/chr{i}/CPMA/{eqtl_file}_meanzscores_chr{i}.gz'
        evalues_file = f'{output}/chr{i}/CPMA/{eqtl_file}_evalues_chr{i}.gz'
        evectors_file = f'{output}/chr{i}/CPMA/{eqtl_file}_Q_chr{i}.gz'
        num_sim = 500000
        sim_file = f'{output}/chr{i}/CPMA/{eqtl_file}_sim{num_sim}_cpma_chr{i}_converted'
        sim_mix_file = f'{output}/chr{i}/CPMA/{eqtl_file}_sim{num_sim}_cpma-mix_chr{i}_converted'
        if not os.path.isfile(sim_file):
            simulate_cmd = f'python ../mixtureModel/simulate_cpma-mix_chunks.py -z {mzscores} -e {evalues_file} -q {evectors_file} -o {sim_file} -m {sim_mix_file} -n {num_sim}'.split(' ')
            simulate = subprocess.Popen(simulate_cmd)
            simulate.wait()
            return_code = simulate.returncode
            chr_error[i] = return_code
    ''' 
    '''
        if len(simulate_cmd) > 0 or (i==22):
            simulate_procs = [ subprocess.Popen(i) for i in simulate_cmd]
            print(simulate_procs)
            for p in simulate_procs:
                p.wait()
                return_code = p.returncode

            simulate_cmd = []       
    
    print('Finished simulation of cpma values')
    '''
    '''
    cpma_file = f'{eqtl_file}_cpma_converted'
   
    num_sim = 500000
    #Compare simulated cpma with observed cpma to get an empirical pvalue for each snp
    compare_cmd = []
    for i in range(1, 23):
        sim_file = f'{output}/chr{i}/CPMA/{eqtl_file}_sim{num_sim}_cpma_chr{i}_converted'
        empirical_file = f'{output}/chr{i}/CPMA/{eqtl_file}_empiricalpvalues_chr{i}_converted'
        compare_cmd.append(f'python calculate_empirical_pvalue.py -s {sim_file} -o {output}/chr{i}/CPMA/{cpma_file}_chr{i} -e {empirical_file}'.split(' '))
        if len(compare_cmd) > 0 or (i==22):
            compare_procs = [ subprocess.Popen(i) for i in compare_cmd]
            print(compare_procs)
            for p in compare_procs:
                p.wait()
            compare_cmd = []
    print('Finished calculating empirical pvalues from cpma for all chr')
    
    #Merge all empirical pvalues of all chr to one file
    all_epvalues = f'{output}/{eqtl_file}_empiricalpvalues_allchr_converted'
    cp_cmd = f'cp {output}/chr1/CPMA/{eqtl_file}_empiricalpvalues_chr1_converted {all_epvalues}'.split(' ')
    cp = subprocess.Popen(cp_cmd).wait()
    for i in range(2, 23):
        append_cmd = f'tail -n +2 -q {output}/chr{i}/CPMA/{eqtl_file}_empiricalpvalues_chr{i}_converted  >> {all_epvalues}'
        subprocess.call(append_cmd, shell=True)
    print('Finished merging all chr empirical pvalues to one file')
    
    for chrm, error in chr_error.items():
        if error != 0:
            print(f'chr{chrm}: Did not finish, Errorcode: {error}')
    '''

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type=str, help='Input folder with expression files and matrix eqtl outputs')
    parser.add_argument('-o', '--output', type=str, help='Destination folder for output')
    input_values = parser.parse_args()

    cpma_pipeline(input_values.input, input_values.output)


if __name__ == '__main__':
    main()
