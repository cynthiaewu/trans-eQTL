import numpy as np
import pandas as pd
import yaml
from scipy.stats import bernoulli
import argparse


def calculate_errors(results, num_targets, num_genes, sig_threshold):
    positive = results[results['p-value'] < sig_threshold]
    #print(positive)
    targets = [f'Gene{i}' for i in range(num_targets)]
    true_pos = len(set(targets) & set(positive['gene']))
    false_neg = num_targets - true_pos
    false_pos = len(positive['gene']) - true_pos
    # print(true_pos, false_neg, false_pos)
    return false_pos/(num_genes-num_targets), false_neg/num_targets  # type I, type II


def get_power(config, input, output):
    with open(config) as f:
        params = yaml.load(f)
        print(params)
    num_genes = params['num_genes']
    num_targets = params['num_targets']
    sig_threshold = params['sig_threshold']
    results = pd.read_csv(input, sep='\t')
    type1, type2 = calculate_errors(results, num_targets, num_genes, sig_threshold)
    power = 1 - type2
    print(f'Power: {power}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=True, help="Input config file with parameters")
    parser.add_argument("-i", "--input", required=True, help="Matrix eqtl results output file as input file")
    parser.add_argument("-o", "--output", required=True, help="Output folder with simulated files")
    params = parser.parse_args()

    get_power(config=params.config,
          input=params.input,
          output=params.output)


if __name__ == "__main__":
    main()
