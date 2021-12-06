import numpy as np
import pandas as pd
import yaml
from scipy.stats import bernoulli
import argparse


def get_pvalues_for_targets(result_file, num_targets):
    results = pd.read_csv(result_file, sep='\t')
    targets = [f'Gene{i}' for i in range(num_targets)]
    results_targets = (results[results.gene.isin(targets)])
    results_targets['gene'] = pd.Categorical(results_targets['gene'], targets)
    #print(results_targets)
    results_targets_sorted = results_targets.sort_values('gene')
    return list(results_targets_sorted['beta']), list(results_targets_sorted['p-value'])


def calculate_power(all_pvalues, sig_threshold):
    pvalues_transpose = np.array(all_pvalues).T
    all_power = []
    for pvalues in pvalues_transpose: # transpose to get pvalues for each beta value used in simulations
        true_pos = len([i for i in pvalues if i < sig_threshold]) 
        total_pos = len(pvalues)
        false_neg = total_pos - true_pos
        type2 = false_neg/total_pos
        all_power.append(1-type2)
    return all_power


def get_power(config, folder, iterations, output):
    with open(config) as f:
        params = yaml.load(f)
        print(params)
    num_genes = params['num_genes']
    num_targets = params['num_targets']
    sig_threshold = params['sig_threshold']
    sim_prefix = 'Simulation'
    all_betas = []
    all_pvalues = []
    for i in range(iterations):
        result_file = f'{folder}/{sim_prefix}_{i}/gene-snp-eqtl.txt'
        betas, pvalues = get_pvalues_for_targets(result_file, num_targets)
        all_betas.append(betas)
        all_pvalues.append(pvalues)
    all_power = calculate_power(all_pvalues, sig_threshold)
    np.savetxt(f'{output}/all_pvalues.txt', all_pvalues)
    np.savetxt(f'{output}/all_betas.txt', all_betas)
    np.savetxt(f'{output}/all_power.txt', all_power)
    print(all_power)
    #type1, type2 = calculate_errors(results, num_targets, num_genes, sig_threshold)
    #power = 1 - type2
    #print(f'Power: {power}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=True, help="Input config file with parameters")
    parser.add_argument("-f", "--folder", required=True, help="Folder with simulation folders which contains matrix eqtl results output file")
    parser.add_argument("-i", "--iterations", type=int, required=True, help="# iterations to simulate genotype and expression files")
    parser.add_argument("-o", "--output", required=True, help="Output folder with simulated files")
    params = parser.parse_args()

    get_power(config=params.config,
          folder=params.folder,
          iterations=params.iterations,
          output=params.output)


if __name__ == "__main__":
    main()
