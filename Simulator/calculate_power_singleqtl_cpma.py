import pandas as pd
import yaml
import argparse


def get_pvalues_for_targets(result_file):
    results = pd.read_csv(result_file, sep='\t')
    return float(results['adj_pvalue'])


def calculate_power(all_pvalues, sig_threshold):
    true_pos = len([i for i in all_pvalues if i < sig_threshold]) 
    total_pos = len(all_pvalues)
    false_neg = total_pos - true_pos
    type2 = false_neg/total_pos
    power = 1-type2
    return power


def get_power(config, folder, iterations):
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
        result_file = f'{folder}/{sim_prefix}_{i}/CPMA/gene-snp-eqtl_cpma_pvalues'
        pvalues = get_pvalues_for_targets(result_file)
        all_pvalues.append(pvalues)
    power = calculate_power(all_pvalues, fdr_sig_threshold)
    calculated = [[num_targets[0], beta_value[0], power]]
    print(calculated)
    power_df = pd.DataFrame(calculated, columns=['#target_genes', 'beta_value', 'power'])
    power_df.to_csv(f'{folder}/power.txt', index=False, sep='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=True, help="Input config file with parameters")
    parser.add_argument("-f", "--folder", required=True, help="Folder with simulation folders which contains simulated data files")
    parser.add_argument("-i", "--iterations", type=int, required=True, help="# iterations to simulate genotype and expression files")
    #parser.add_argument("-o", "--output", required=True, help="Output folder to write power analysis files")
    params = parser.parse_args()

    get_power(config=params.config,
          folder=params.folder,
          iterations=params.iterations)
         # output=params.output)


if __name__ == "__main__":
    main()
