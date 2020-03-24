import pandas as pd
import numpy as np
import argparse


def filter_genes(input, output):
    #data = pd.read_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr2/Clean_expression_chr2_inter_filtered.tsv', sep='\t')
    data = pd.read_csv(input, sep='\t')
    data = data.set_index('Unnamed: 0')
    indices = np.where(np.sum(data > 0.1, axis=1) >= 10)[0]
    filtered_data = data.iloc[indices,:]
    filtered_data.to_csv(output, header=True, sep='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input expression file with filtered protein coding genes")
    parser.add_argument("-o", "--output", required=True, help="Ouptput file with further filtered genes, rpkm > 0.1 in 10 samples")
    params = parser.parse_args()
    filter_genes(params.input, params.output)


if __name__ == "__main__":
    main()
