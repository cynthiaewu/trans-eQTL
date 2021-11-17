import pandas as pd
import numpy as np
from numpy import genfromtxt
import math
import argparse


def calculate_cpma(input, num_genes, num_snps, output):
    #pvalues = genfromtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_pvalues_nofdr.csv', delimiter='\t', skip_header=1)
    pvalues = genfromtxt(input, delimiter='\t', skip_header=1)
    #np.negative(np.log(pvalues[1:]))
    #num_snps = len(pvalues)
    #num_genes = len(pvalues[1])-1
    if num_snps == 1:
        pvalues = pvalues.reshape(1, -1)
    cpma = []
    for i in range(num_snps):
        likelihood = np.mean(np.negative(np.log(pvalues[i][1:])))
        value = -2 * ((((likelihood - 1) * num_genes)/likelihood) - num_genes*np.log(likelihood))
        if math.isnan(value):
            value = 0
        cpma.append(value)

    snps = pd.read_csv(input, usecols=[0], sep='\t', names = ['snp'], skiprows = 1)
    snps.insert(column = 'cpma', loc = 1, value = cpma)
    #snps.to_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_cpma_nofdr.csv', index=False, header=True, sep='\t')
    snps.to_csv(output, index=False, header=True, sep='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input pvalues no fdr file")
    parser.add_argument("-g", "--genes", type=int, required=True, help="# genes")
    parser.add_argument("-s", "--snps", type=int, required=True, help="# snps")
    parser.add_argument("-o", "--output", required=True, help="Ouptput file with cpma values")
    params = parser.parse_args()
    calculate_cpma(params.input, params.genes, params.snps, params.output)


if __name__ == "__main__":
    main()
