import pandas as pd
import numpy as np
import scipy.optimize, scipy.stats
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
    top_numgenes = math.ceil(x*num_genes)
    cpma = []
    for i in range(num_snps):
        all_values = pvalues[i][1:]
        all_values.sort()
        all_values = -np.log(all_values)
        top_values = all_values[:top_numgenes]
        bordernum = top_values[-1]
        #print(top_values)
        def log_likelihood_neg(x):
            #print(x)
            log_likelihood = top_numgenes*np.log(x) - \
                             x*np.sum(top_values) + \
                             (num_genes-top_numgenes)*np.log(1-np.exp(-x*bordernum))
            return -log_likelihood
        if math.isclose(x, 1.0):
            max_likelihood = 1/np.mean(all_values)
        else:
            max_likelihood = scipy.optimize.minimize(log_likelihood_neg,
                                                     x0=1,
                                                     bounds=((10**(-10), None),))
            max_likelihood = max_likelihood.x[0]
        #print(max_likelihood)
        term1 = (1-max_likelihood)*np.sum(top_values)
        term2 = top_numgenes*np.log(max_likelihood)
        term3 = (num_genes-top_numgenes)*np.log(1-math.exp(-max_likelihood*bordernum))
        term4 = -(num_genes-top_numgenes)*np.log(1-math.exp(-bordernum))
        #print(term1, term2, term3, term4)
        cpmax_value = 2*(term1 + term2 + term3 + term4)
        pvalue = 1 - scipy.stats.chi2.cdf(cpmax_value, 1)
        #print(cpmax_value, pvalue)
        if math.isnan(cpmax_value):
            cpmax_value = 0
        cpma.append(cpmax_value)

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
