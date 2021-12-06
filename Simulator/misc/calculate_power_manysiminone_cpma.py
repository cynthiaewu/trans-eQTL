import numpy as np
import pandas as pd
import yaml
from scipy.stats import bernoulli
import argparse


def get_pvalues_for_targets(results, count, rep):
    targets = [f'SNP{i}' for i in range(count, count+rep)]
    results_targets = (results[results.snp.isin(targets)])
    return list(results_targets['adj_pvalue'])


def calculate_power(pvalues, sig_threshold):
    true_pos = len([i for i in pvalues if i < sig_threshold]) 
    total_pos = len(pvalues)
    false_neg = total_pos - true_pos
    type2 = false_neg/total_pos
    return (1-type2)


def get_power(config, result_file, output):
    with open(config) as f:
        params = yaml.load(f)
        print(params)
    num_targets = params['num_targets']
    fixed_betas = params['beta_value']
    rep = params['rep']
    sig_threshold = params['sig_threshold']
    #fdr_sig_threshold = sig_threshold/rep
    results = pd.read_csv(result_file, sep='\t')
    # counter to keep track of the ith snp already evaluated
    # snp numbering is based for all #target, for all beta value, #snps for each = #rep
    count = 0
    calculated = []
    for num_t in num_targets:
        for beta_value in fixed_betas:
            pvalues = get_pvalues_for_targets(results, count, rep)
            power = calculate_power(pvalues, sig_threshold)
            calculated.append([num_t, beta_value, power])
            count += rep
    power_df = pd.DataFrame(calculated, columns=['#target_genes', 'beta_value', 'power'])
    power_df.to_csv(output, index=False, sep='\t')
    #print(power_df)
    #type1, type2 = calculate_errors(results, num_targets, num_genes, sig_threshold)
    #power = 1 - type2
    #print(f'Power: {power}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=True, help="Input config file with parameters")
    parser.add_argument("-i", "--result_file", required=True, help="File with cpma, pvalue, adj pvalue of the snps")
    parser.add_argument("-o", "--output", required=True, help="Output file to write power analysis file")
    params = parser.parse_args()

    get_power(config=params.config,
          result_file=params.result_file,
          output=params.output)


if __name__ == "__main__":
    main()
