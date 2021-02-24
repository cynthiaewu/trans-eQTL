import pandas as pd
import numpy as np
import yaml
import argparse

def calculate_type1(nullpvals, sig_threshold):
    false_pos = len([i for i in nullpvals if i < sig_threshold])
    total_neg = len(nullpvals)
    #print(false_pos, total_neg)
    type1 = false_pos/total_neg
    return type1


def iterate_files(config, folder, iterations, num_snps_real):
    with open(config) as f:
        params = yaml.load(f)
        print(params)
    num_snps = params['num_snps']
    null_snps = params['num_nullsnps']
    num_targets = params['num_targets']
    beta_value = params['beta_value']
    sig_threshold = params['sig_threshold']
    fdr_sig_threshold = sig_threshold/num_snps_real
    eqtl = ['SNP' + str(i) for i in range(num_snps-null_snps)]
    type1= []
    for i in range(iterations):
        #result_file = f'{folder}/Simulation_{i}/CPMAx_PEER/gene-snp-eqtl_cpmax_pvalues_1.0'
        result_file = f'{folder}/Simulation_{i}/CPMAx/gene-snp-eqtl_cpmax_pvalues_1.0'
        #result_file = f'{folder}/Simulation_{i}/CPMA/gene-snp-eqtl_empiricalpvalues_topx_1.0'
        #result_file = f'{folder}/Simulation_{i}/expressionPCs/gene-snp-eqtl_PCs_cpmax_pvalues_1.0'
        results = pd.read_csv(result_file, sep='\t')
        nullpvals = np.array(results.loc[~results['snp'].isin(eqtl)]['pvalue'])
        type1.append(calculate_type1(nullpvals, fdr_sig_threshold))
    #print(type1)
    np.savetxt(f'{folder}/type1_error_rate_{num_snps_real}tests_uncorrectedPEER.txt', type1, delimiter='\t')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=True, help="Input config file with parameters")
    parser.add_argument("-f", "--folder", required=True, help="Folder with simulation folders which contains simulated data files")
    parser.add_argument("-i", "--iterations", type=int, required=True, help="# iterations to simulate genotype and expression files")
    parser.add_argument("-n", "--num_snps_real", default=10000, type=int, help="# snps to use for multiple hypothesis correction (based on real data)")
    params = parser.parse_args()
    iterate_files(params.config, params.folder, params.iterations, params.num_snps_real)

if __name__ == '__main__':
    main()
