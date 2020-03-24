import pandas as pd
import numpy as np
import argparse


def get_european(input, european, dtype, output):
    #data = pd.read_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr2/Clean_expression_chr2_inter_filtered.tsv', sep='\t')
    with open(european) as f:
        euro_samples = f.read().splitlines()
    
    data = pd.read_csv(input, sep='\t')
    if dtype == 0:
        data = data.set_index('chrom_start')
    if dtype == 1:
        data = data.set_index('Unnamed: 0')
    euro_data = data[euro_samples]
    euro_data.to_csv(output, sep='\t')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input data file")
    parser.add_argument("-e", "--european", required=True, help="Input european samples file")
    parser.add_argument("-t", "--dtype", required=True, type = int, help="Data type, 0 if genotype file, 1 if expression file")
    parser.add_argument("-o", "--output", required=True, help="Ouptput data file with only European samples")
    params = parser.parse_args()
    get_european(params.input, params.european, params.dtype, params.output)


if __name__ == "__main__":
    main()
