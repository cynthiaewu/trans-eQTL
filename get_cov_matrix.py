import pandas as pd
import numpy as np
import argparse


def getCov(input, output):
    data = pd.read_csv(input, sep='\t', index_col=0)
    genes = list(data.columns)
    scores = np.array(data)
    cov = np.cov(np.transpose(scores))
    cov_matrix = pd.DataFrame(cov, columns=genes)
    #cov_matrix.to_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_cov.csv', index=True, header=True, sep='\t')
    cov_matrix.to_csv(output, index=False, header=True, sep='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input file with zscores from matrix eQTL")
    parser.add_argument("-o", "--output", required=True, help="Ouptput file with simulated zscores")
    params = parser.parse_args()
    getCov(params.input, params.output)


if __name__ == "__main__":
    main()
