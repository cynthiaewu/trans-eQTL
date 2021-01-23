import pandas as pd
import numpy as np
import argparse

def getValues(input, p_out, z_out): 
    #data = pd.read_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls.csv', sep='\t')
    data = pd.read_csv(input, sep='\t')
    zscores = data.pivot_table(index='SNP', columns='gene', values='t-stat')
    pvalues = data.pivot_table(index='SNP', columns='gene', values='p-value')

    pvalues.to_csv(p_out, index=True, header=True, sep='\t')
    zscores.to_csv(z_out, index=True, header=True, sep='\t')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input matrix eQTL file")
    parser.add_argument("-p", "--pvalues_out", required=True, help="Output file with matrix of pvalues for every snp and gene")
    parser.add_argument("-z", "--zscores_out", required=True, help="Output file with matrix of zscores for every snp and gene")
    params = parser.parse_args()
    getValues(params.input, params.pvalues_out, params.zscores_out)


if __name__ == "__main__":
    main()
