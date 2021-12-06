import pandas as pd
import numpy as np
import scipy.optimize, scipy.stats
from numpy import genfromtxt
import scipy.stats as stats 
import math
import argparse


def likelihood_ratio_test(pvals):
    fit_alpha, fit_loc, fit_beta=stats.gamma.fit(pvals, floc=0)
    alt_lklh = stats.gamma.logpdf(pvals, fit_alpha, fit_loc, fit_beta).sum()
    null_lklh = stats.gamma.logpdf(pvals, 1,0,1).sum()

    test_stat = 2*(-null_lklh + alt_lklh)
    pvalue = 1 - stats.chi2.cdf(test_stat, 2)
    return test_stat, pvalue


def get_pvalue(input, output):
    #pvalues = genfromtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_pvalues_nofdr.csv', delimiter='\t', skip_header=1)
    pvalues = genfromtxt(input, delimiter='\t', skip_header=1)
    #np.negative(np.log(pvalues[1:]))
    if len(pvalues.shape) == 1:
      pvalues = pvalues.reshape(1, -1)
    #print(len(pvalues.shape))
    num_snps = len(pvalues)
    #num_genes = len(pvalues[1])-1
    num_genes = len(pvalues[0])-1
    all_metapvals = []
    all_teststat = []
    for i in range(num_snps):
        all_values = pvalues[i][1:]
        all_values = -np.log(all_values)
        test_stat, meta_pval = likelihood_ratio_test(all_values)
        all_teststat.append(test_stat)
        all_metapvals.append(meta_pval)
    snps = pd.read_csv(input, usecols=[0], sep='\t', names = ['snp'], skiprows = 1)
    snps.insert(column = 'test_stat', loc = 1, value = all_teststat)
    snps.insert(column = 'gamma_pval', loc = 2, value = all_metapvals)
    #snps.to_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_cpma_nofdr.csv', index=False, header=True, sep='\t')
    snps.to_csv(output, index=False, header=True, sep='\t')
    


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input pvalues no fdr file")
    parser.add_argument("-o", "--output", required=True, help="Ouptput file with cpma values")
    params = parser.parse_args()
    get_pvalue(params.input, params.output)


if __name__ == "__main__":
    main()
