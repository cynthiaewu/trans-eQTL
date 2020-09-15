import pandas as pd
import numpy as np
from numpy import genfromtxt
import math
import argparse


def calculate_cpma_topx(input, x, output):
    #pvalues = genfromtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_pvalues_nofdr.csv', delimiter='\t', skip_header=1)
    pvalues = genfromtxt(input, delimiter='\t', skip_header=1)
    #np.negative(np.log(pvalues[1:]))
    if len(pvalues.shape) == 1:
      pvalues = pvalues.reshape(1, -1)
    #print(len(pvalues.shape))
    num_snps = len(pvalues)
    #num_genes = len(pvalues[1])-1
    num_genes = len(pvalues[0])-1
    top_ngenes = math.ceil(x*num_genes)
    cpma = []
    for i in range(num_snps):
        top_values = sorted(pvalues[i][1:])[:top_ngenes]
        #likelihood = np.mean(np.negative(np.log(pvalues[i][1:])))
        likelihood = np.mean(np.negative(np.log(top_values)))
        value = -2 * ((((likelihood - 1) * num_genes)/likelihood) - num_genes*np.log(likelihood))
        if math.isnan(value):
            value = 0
        cpma.append(value)

    #print(cpma)
    snps = pd.read_csv(input, usecols=[0], sep='\t', names = ['snp'], skiprows = 1)
    snps.insert(column = 'cpma', loc = 1, value = cpma)
    #snps.to_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_cpma_nofdr.csv', index=False, header=True, sep='\t')
    snps.to_csv(output, index=False, header=True, sep='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input pvalues no fdr file")
    parser.add_argument("-x", "--xpercent", default=0.1, type=float, help="Top x percent of genes to be used for cpma calculation")
    parser.add_argument("-o", "--output", required=True, help="Ouptput file with cpma values")
    params = parser.parse_args()
    calculate_cpma_topx(params.input, params.xpercent, params.output)


if __name__ == "__main__":
    main()
