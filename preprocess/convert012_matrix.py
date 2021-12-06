import argparse
import pandas as pd
import numpy as np


def convert_matrix(fname, output):
    #data = np.loadtxt('/storage/cynthiawu/trans_eQTL/GTex.v8/Run4_080221/chr1/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv_filter3_snplist_chr1.012', delimiter='\t')
    data = np.loadtxt(fname, delimiter='\t')
    data = np.delete(data, 0, axis=1)

    indv = []
    with open(f'{fname}.indv', 'r') as f:
        for line in f:
            indv.append(line.strip())

    pos = []
    with open(f'{fname}.pos', 'r') as f:
        for line in f:
            pos.append(line.strip().replace('\t', '_'))
    genotype = pd.DataFrame(data.T, index=pos, columns=indv)
    genotype.to_csv(output, sep='\t')

        


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fname", required=True, help="Input file with .012 output from vcftools")
    parser.add_argument("-o", "--output", required=True, help="Output file with reformat genotype matrix")
    params = parser.parse_args()
    convert_matrix(params.fname, params.output)


if __name__ == "__main__":
    main()
