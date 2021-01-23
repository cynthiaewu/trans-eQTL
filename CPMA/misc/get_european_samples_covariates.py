import pandas as pd
import numpy as np
import argparse


def get_european(input, european, output):
    #data = pd.read_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr2/Clean_expression_chr2_inter_filtered.tsv', sep='\t')
    with open(european) as f:
        euro_samples = f.read().splitlines()
    euro_samples_format = []
    for s in euro_samples:
        euro_samples_format.append(s.replace('.', '-'))

    #cov_data = pd.read_csv('/storage/cynthiawu/trans_eQTL/gtex_covariates.txt', sep='\t', )
    cov_data = pd.read_csv(input, sep='\t', )
    cov_data.set_index('SUBJID',inplace=True)
    cov_dataT = cov_data.T
    euro_data = cov_dataT[euro_samples_format]
    euro_data = euro_data.reset_index()
    euro_data = euro_data.rename(columns={'index': 'SUBJID'})
    #euro_data.to_csv('/storage/cynthiawu/trans_eQTL/European_samples_covariates.txt', sep='\t', index = False)
    euro_data.to_csv(output, sep='\t', index = False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input data file")
    parser.add_argument("-e", "--european", required=True, help="Input european samples file")
    parser.add_argument("-o", "--output", required=True, help="Ouptput data file with only European samples")
    params = parser.parse_args()
    get_european(params.input, params.european, params.output)


if __name__ == "__main__":
    main()
