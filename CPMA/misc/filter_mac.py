import pandas as pd
import numpy as np
from collections import Counter
import argparse


def filter_genes(input, output):
    data = pd.read_csv(input, sep='\t')
    discard = []
    for index, row in data.iterrows():
        genotypes = list(row[1:])
        # get the minor allele count by taking top 2 of the counts of unique genotypes, least
        mac = Counter(genotypes).most_common(2)[-1][1]
        if mac < 2:
            discard.append(index)
    data_mac = data.drop(data.index[discard])

    data_mac.to_csv(output, header=True, index=False, sep='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input expression file")
    parser.add_argument("-o", "--output", required=True, help="Ouptput filtered genotype file by mac (after intersecting samples)")
    params = parser.parse_args()
    filter_genes(params.input, params.output)


if __name__ == "__main__":
    main()
