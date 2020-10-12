import pandas as pd
import numpy as np
from scipy import stats
import yaml
import argparse


def get_pvalues_for_targets(result_file, method):
    results = pd.read_csv(result_file, sep='\t')
    if method==3 or method==4:
        samplesize = len(results['p-value'])
        pvalues = results['p-value']
    if method==3:
        return min(float(pvalues[0])*samplesize, 1)
    if method==4:
        return  stats.kstest((list(pvalues)), 'uniform')[1]
    else:
        return float(results['adj_pvalue'])


def calculate_power(all_pvalues, sig_threshold):
    true_pos = len([i for i in all_pvalues if i < sig_threshold]) 
    total_pos = len(all_pvalues)
    false_neg = total_pos - true_pos
    type2 = false_neg/total_pos
    power = 1-type2
    return power


def get_power(config, method, topx, folder, iterations):
    with open(config) as f:
        params = yaml.load(f)
        print(params)
    num_snps = params['num_snps']
    num_targets = params['num_targets']
    beta_value = params['beta_value']
    sig_threshold = params['sig_threshold']
    fdr_sig_threshold = sig_threshold/num_snps
    sim_prefix = 'Simulation'
    all_pvalues = []
    for i in range(iterations):
        if method==0:
            result_file = f'{folder}/{sim_prefix}_{i}/CPMA/gene-snp-eqtl_cpma_pvalues_fixed'
        if method==1:
            result_file = f'{folder}/{sim_prefix}_{i}/CPMAx/gene-snp-eqtl_cpmax_pvalues_{topx}'
        # cpmax on PC expression matrix 
        if method==2:
            result_file = f'{folder}/{sim_prefix}_{i}/expressionPCs/gene-snp-eqtl_PCs_cpmax_pvalues_{topx}'
        # only matrix eqtl results on PC expression matrix
        if method==3 or method==4:
            result_file = f'{folder}/{sim_prefix}_{i}/expressionPCs/gene-snp-eqtl_PCs'
        pvalues = get_pvalues_for_targets(result_file, method)
        all_pvalues.append(pvalues)
    power = calculate_power(all_pvalues, fdr_sig_threshold)
    calculated = [[num_targets[0], beta_value[0], power]]
    print(f'calculated, {folder}')
    power_df = pd.DataFrame(calculated, columns=['#target_genes', 'beta_value', 'power'])
    if method==0:
        power_df.to_csv(f'{folder}/power_fixed.txt', index=False, sep='\t')
    if method==1:
        power_df.to_csv(f'{folder}/power_cpmax_{topx}.txt', index=False, sep='\t')
    if method==2:
        power_df.to_csv(f'{folder}/power_PCs_cpmax_{topx}.txt', index=False, sep='\t')
    if method==3:
        power_df.to_csv(f'{folder}/power_PCs.txt', index=False, sep='\t')
    if method==4:
        power_df.to_csv(f'{folder}/power_PCs_kstest.txt', index=False, sep='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=True, help="Input config file with parameters")
    parser.add_argument("-m", "--method", type=int, required=True, help="0 for cpma, 1 for cpma_topx, or 2 for cpma_topx_PCs")
    parser.add_argument("-x", "--topx", default=0.1, type=float, help="Top x percent of genes to be used for cpma calculation")
    parser.add_argument("-f", "--folder", required=True, help="Folder with simulation folders which contains simulated data files")
    parser.add_argument("-i", "--iterations", type=int, required=True, help="# iterations to simulate genotype and expression files")
    #parser.add_argument("-o", "--output", required=True, help="Output folder to write power analysis files")
    params = parser.parse_args()

    get_power(config=params.config,
          method=params.method,
          topx=params.topx,
          folder=params.folder,
          iterations=params.iterations)
         # output=params.output)


if __name__ == "__main__":
    main()
