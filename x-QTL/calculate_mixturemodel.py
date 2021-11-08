import pandas as pd
import numpy as np
import scipy.optimize, scipy.stats
from numpy import genfromtxt
import math
import argparse


def log_likelihood_neg(t, L, pvals):
    pvals = np.array(pvals)
    return -np.sum(np.log((1 - t) * np.exp(-pvals) + t * 1/L * np.exp(-1/L*pvals)))


def likelihood_ratio_test(pvals, trueT=None, trueL=None):
    assert (trueT is None) == (trueL is None)
    null_lklh = log_likelihood_neg(0, 1, pvals)
    if trueT is None:
        results = scipy.optimize.minimize(lambda tL: log_likelihood_neg(*tL, pvals=pvals),
                                           method='L-BFGS-B',
                                           x0=(0.5, 1),
                                           bounds=((10**(-5), 1),(10**(-5), None)))
        alt_lklh = results['fun']
        x = results['x']
        bestT, bestL = x[0], x[1]
    else:
        alt_lklh = log_likelihood_neg(t=trueT, L=trueL, pvals=pvals)
    test_stat = -2*(-null_lklh + alt_lklh)
    
    #pvalue = 1 - scipy.stats.chi2.cdf(test_stat, 1)
    return test_stat, bestT, bestL


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
    #all_metapvals = []
    all_teststats = []
    all_T = []
    all_L = []
    for i in range(num_snps):
        all_values = pvalues[i][1:]
        all_values = np.where(all_values > 10**(-150), all_values, 10**(-150))
        all_values = -np.log(all_values)
        test_stat, best_T, best_L = likelihood_ratio_test(all_values)
        #all_metapvals.append(meta_pval)
        all_teststats.append(test_stat)
        all_T.append(best_T)
        all_L.append(best_L)
    snps = pd.read_csv(input, usecols=[0], sep='\t', names = ['snp'], skiprows = 1)
    #snps.insert(column = 'mixture_pval', loc = 1, value = all_metapvals)
    snps.insert(column = 'mixture_tstat', loc = 1, value = all_teststats)
    snps.insert(column = 'predicted_T', loc = 2, value = all_T)
    snps.insert(column = 'predicted_L', loc = 3, value = all_L)
    #snps.to_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_cpma_nofdr.csv', index=False, header=True, sep='\t')
    snps.to_csv(output, index=False, header=True, sep='\t')
    

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input pvalues no fdr file")
    parser.add_argument("-o", "--output", required=True, help="Ouptput file with cpma-mix values")
    params = parser.parse_args()
    get_pvalue(params.input, params.output)


if __name__ == "__main__":
    main()
