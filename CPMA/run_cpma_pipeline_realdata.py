import sys
import subprocess
import os
import argparse

def cpma_pipeline(input_folder, output):
    eqtl_file = f'gene-snp_eqtls'
    
    if not os.path.isdir(output):
        os.mkdir(output)
    
    #for i in range(1, 23):
    for i in [1]:
        chr_folder = os.path.join(output, f'chr{i}')
        if not os.path.isdir(chr_folder):
            os.mkdir(chr_folder)
    
    # Obtain zscores and pvalues in a snp by gene matrix format from matrix eQTL output
    #for i in range(1, 23):
    for i in [1]:
        pvalue_file = f'{eqtl_file}_pvalue'
        zscore_file = f'{eqtl_file}_zscore'
    
        values_cmd = f'python get_values.py -i {input_folder}/chr{i}_matrixeQTL_raw.tsv.gz -z {output}/chr{i}/{zscore_file}_chr{i}'.split(' ')
        values = subprocess.Popen(values_cmd).wait()
        convert_cmd = f'python tstat_to_pvalue.py -i {output}/chr{i}/{zscore_file}_chr{i} -p {output}/chr{i}/{pvalue_file}_chr{i}'.split(' ')
        convert = subprocess.Popen(convert_cmd).wait()
        print(f'Finished getting pvalues and zscores for chr{i}')
    print('Finished getting pvalues and zscores for all chr')

    #Calculate cpma values for all chr with matrix eqtl pvalues output
    cpma_file = f'{eqtl_file}_cpma'
    
    #for i in range(1, 23):
    for i in [1]:
        cpma_cmd = f'python calculate_cpma.py -i {output}/chr{i}/{pvalue_file}_chr{i} -o {output}/chr{i}/{cpma_file}_chr{i}'.split(' ')
        cpma = subprocess.Popen(cpma_cmd).wait()
        print(f'Finished calculating cpma for chr{i}')
    print('Finished calculating cpma for all chr')
    
    '''
    #Combine zscores of all chr in one file
    all_out = f'{output}/{eqtl_file}_allchr'
    all_zscores = f'{all_out}_zscores'
    
    cp_cmd = f'cp {output}/chr1/{zscore_file}_chr1 {all_zscores}'.split(' ')
    cp = subprocess.Popen(cp_cmd).wait()
    for i in range(2, 23):
        append_cmd = f'tail -n +2 -q {output}/chr{i}/{zscore_file}_chr{i}  >> {all_zscores}'
        subprocess.call(append_cmd, shell=True)
    print('Finished merging all chr zscores to one file')
    '''
    
    cov_cmd = []
    #for i in range(1, 23):
    for i in [1]:
        #cov_matrix = f'{output}/chr{i}/{eqtl_file}_cov_chr{i}'
        mzscores = f'{output}/chr{i}/{eqtl_file}_meanzscores_chr{i}.gz'
        evalues_file = f'{output}/chr{i}/{eqtl_file}_evalues_chr{i}.gz'
        evectors_file = f'{output}/chr{i}/{eqtl_file}_Q_chr{i}.gz'
        cov_cmd.append(f'python calculate_cov_meanzscores_edecomposition.py -i {output}/chr{i}/{eqtl_file}_zscore_chr{i} -m {mzscores} -e {evalues_file} -q {evectors_file}'.split(' '))

        if len(cov_cmd) > 1 or (i==22):
            cov_procs = [ subprocess.Popen(i) for i in cov_cmd]
            print(cov_procs)
            for p in cov_procs:
                p.wait()
            cov_cmd = []
        #cov = subprocess.Popen(cov_cmd).wait()
    print('Finished calculating mean zscores and eigendecomposition')
    
    '''
    #Calculate the gene covariance matrix and mean zscores for genes
    cov_cmd = []
    for i in range(1, 23):
        cov_matrix = f'{output}/chr{i}/{eqtl_file}_cov_chr{i}'
        mzscores = f'{output}/chr{i}/{eqtl_file}_meanzscores_chr{i}'
        cov_cmd.append(f'python calculate_cov_meanzscores.py -i {output}/chr{i}/{eqtl_file}_zscore_chr{i} -c {cov_matrix} -m {mzscores}'.split(' '))
        if len(cov_cmd) > 3 or (i==22):
            cov_procs = [ subprocess.Popen(i) for i in cov_cmd]
            print(cov_procs)
            for p in cov_procs:
                p.wait()
            cov_cmd = []
        #cov = subprocess.Popen(cov_cmd).wait()
    print('Finished calculating cov matrix and mean zscores')
    
    
    #Calculate the eigendecomposition
    eigen_cmd = []
    for i in range(1, 23):
        cov_matrix = f'{output}/chr{i}/{eqtl_file}_cov_chr{i}'
        evalues_file = f'{output}/chr{i}/{eqtl_file}_evalues_chr{i}'
        evectors_file = f'{output}/chr{i}/{eqtl_file}_Q_chr{i}'
        eigen_cmd.append(f'python calculate_edecomposition.py -c {cov_matrix} -e {evalues_file} -q {evectors_file}'.split(' '))
        if len(eigen_cmd) > 0 or (i==22):
            eigen_procs = [ subprocess.Popen(i) for i in eigen_cmd]
            print(eigen_procs)
            for p in eigen_procs:
                p.wait()
            eigen_cmd = []
        #eigen = subprocess.Popen(eigen_cmd).wait()
    print('Finished calculating eigendecomposition')
    '''
    
    #Simulate cpma values from normal distribution with gene covariance matrix and mean zscores
    simulate_cmd = []
    chr_error = {}
    #for i in range(1, 23):
    for i in [1]:
    #for i in [10, 11, 12, 16, 17, 21]:
        mzscores = f'{output}/chr{i}/{eqtl_file}_meanzscores_chr{i}.gz'
        evalues_file = f'{output}/chr{i}/{eqtl_file}_evalues_chr{i}.gz'
        evectors_file = f'{output}/chr{i}/{eqtl_file}_Q_chr{i}.gz'
        num_sim = 500000
        sim_file = f'{output}/chr{i}/{eqtl_file}_sim{num_sim}_cpma_chr{i}'
        simulate_cmd = f'python simulate_cpma_chunks.py -z {mzscores} -e {evalues_file} -q {evectors_file} -o {sim_file} -n {num_sim}'.split(' ')
        simulate = subprocess.Popen(simulate_cmd)
        simulate.wait()
        return_code = simulate.returncode
        chr_error[i] = return_code
    
    
    
    '''
        if len(simulate_cmd) > 0 or (i==22):
            simulate_procs = [ subprocess.Popen(i) for i in simulate_cmd]
            print(simulate_procs)
            for p in simulate_procs:
                p.wait()
                return_code = p.returncode

            simulate_cmd = []
    '''        
    
    
    print('Finished simulation of cpma values')
    
    num_sim = 500000
    #Compare simulated cpma with observed cpma to get an empirical pvalue for each snp
    compare_cmd = []
    #for i in range(1, 23):
    for i in [1]:
        sim_file = f'{output}/chr{i}/{eqtl_file}_sim{num_sim}_cpma_chr{i}'
        empirical_file = f'{output}/chr{i}/{eqtl_file}_empiricalpvalues_chr{i}'
        compare_cmd.append(f'python calculate_empirical_pvalue.py -s {sim_file} -o {output}/chr{i}/{cpma_file}_chr{i} -e {empirical_file}'.split(' '))
        if len(compare_cmd) > 0 or (i==22):
            compare_procs = [ subprocess.Popen(i) for i in compare_cmd]
            print(compare_procs)
            for p in compare_procs:
                p.wait()
            compare_cmd = []
    print('Finished calculating empirical pvalues from cpma for all chr')
    '''
    #Merge all empirical pvalues of all chr to one file
    all_epvalues = f'{output}/{eqtl_file}_empiricalpvalues_allchr'
    cp_cmd = f'cp {output}/chr1/{eqtl_file}_empiricalpvalues_chr1 {all_epvalues}'.split(' ')
    cp = subprocess.Popen(cp_cmd).wait()
    for i in range(2, 23):
        append_cmd = f'tail -n +2 -q {output}/chr{i}/{eqtl_file}_empiricalpvalues_chr{i}  >> {all_epvalues}'
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
