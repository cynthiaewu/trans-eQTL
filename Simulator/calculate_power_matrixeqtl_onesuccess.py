import pandas as pd
import numpy as np
import yaml
import argparse


def get_pvalues_for_targets(result_file, num_targets):
    results = pd.read_csv(result_file, sep='\t')
#     print(len([i for i in list(results['p-value']) if i < 0.05]))

    targets = [f'Gene{i}' for i in range(num_targets)]
    results_targets = (results[results.gene.isin(targets)])
    pvalues = list(results_targets['p-value'])
    # get the smallest pvalue, this will be the snp-gene pair used to measure "success" of detecting if snp is trans-eqtl
    return min(pvalues)


def calculate_power(pvalues, sig_threshold):
    true_pos = len([i for i in pvalues if i < sig_threshold])
    total_pos = len(pvalues)
    false_neg = total_pos - true_pos
    type2 = false_neg/total_pos
    power = 1-type2
    return power


def get_power(config, folder, num_snps_real, iterations):
    with open(config) as f:
        params = yaml.load(f)
        print(params)
    num_genes = params['num_genes']
    num_targets = params['num_targets']
    beta_value = params['beta_value']
    sig_threshold = params['sig_threshold']
    # iterations = num_snps
    correction_threshold = sig_threshold/(num_genes*num_snps_real)
    sim_prefix = 'Simulation'
    pvalues = []
    for i in range(iterations):
        result_file = f'{folder}/{sim_prefix}_{i}/CPMA/gene-snp-eqtl'
        pvalues.append(get_pvalues_for_targets(result_file, num_targets[0]))
    power = calculate_power(pvalues, correction_threshold)
    calculated = [[num_targets[0], beta_value[0], power]]
    power_df = pd.DataFrame(calculated, columns=['#target_genes', 'beta_value', 'power_adjusted'])
    power_df.to_csv(f'{folder}/power_matrix-eqtl_onesuccess_realnumsnps.txt', index=False, sep='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=True, help="Input config file with parameters")
    parser.add_argument("-f", "--folder", required=True, help="Folder with simulation folders which contains matrix eqtl results output file")
    parser.add_argument("-n", "--num_snps_real", default=10000, type=int, help="# snps to use for multiple hypothesis correction (based on real data)")
    parser.add_argument("-i", "--iterations", type=int, required=True, help="# iterations to simulate genotype and expression files")
    params = parser.parse_args()

    get_power(config=params.config,
          folder=params.folder,
          num_snps_real=params.num_snps_real,
          iterations=params.iterations)


if __name__ == "__main__":
    main()

