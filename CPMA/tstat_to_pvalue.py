import pandas as pd
import numpy as np
from scipy import stats
import argparse

def getValues(input, p_out): 
    zscores = pd.read_csv(input, sep='\t', index_col='SNP')
    pvalues = 2*stats.norm.cdf(-np.abs(zscores))
    pvalues_df = pd.DataFrame(pvalues, index=zscores.index, columns=zscores.columns)

    pvalues_df.to_csv(p_out, index=True, header=True, sep='\t')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_zscores", required=True, help="Input matrix of zscores for every snp and gene from matrix eqtl")
    parser.add_argument("-p", "--pvalues_out", required=True, help="Output file with matrix of pvalues for every snp and gene")
    params = parser.parse_args()
    getValues(params.input_zscores, params.pvalues_out)


if __name__ == "__main__":
    main()
