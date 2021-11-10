import pandas as pd
import numpy as np
import scipy.optimize, scipy.stats
from numpy import genfromtxt
from scipy.stats import gamma as gamma_func
from scipy.special import gamma, factorial
import math
import argparse


def log_likelihood_neg_gamma(t, k, theta, pvals):
    gamma_pdf_scipy = gamma_func.pdf(pvals, k, 0, theta)
    #gamma_pdf = (pvals**(k-1)*np.exp(-pvals/theta))/((theta**k)*gamma(k))
    return -np.sum(np.log((1-t) * np.exp(-pvals) + t * gamma_pdf_scipy))

'''
def fit_gamma(pvals):
    fit_alpha, fit_loc, fit_beta=stats.gamma.fit(pvals, floc=0)
    alt_lklh = stats.gamma.logpdf(pvals, fit_alpha, fit_loc, fit_beta).sum()
    null_lklh = stats.gamma.logpdf(pvals, 1,0,1).sum()


def log_likelihood_neg_gamma(t, alpha, beta, pvals):
    pvals = np.array(pvals)
    n = len(pvals)
    s1 = np.log(pvals).sum()
    s2 = np.sum(pvals)
    gamma = n*(alpha*np.log(beta) - np.log(gamma(a))) + (alpha-1)*s1 - beta*s2

    #s2 = sc.log1p(-data).sum()
    return -np.sum(np.log((1 - t) * np.exp(-pvals) + t * 1/L * np.exp(-1/L*pvals)))
'''


def likelihood_ratio_test(pvals):
    null_lklh = log_likelihood_neg_gamma(0, 1, 1, pvals)
    #print(null_lklh)
    alt_lklh = scipy.optimize.minimize(lambda tktheta: log_likelihood_neg_gamma(*tktheta, pvals=pvals),
                                       method='L-BFGS-B',
                                       x0=(0.5, 1, 1),
                                       bounds=(
                                           (10**(-5), 1-10**(-5)),
                                           (10**(-5), None),
                                           (10**(-5), None)
                                       )
                                       )['fun']
    #print(alt_lklh)
    test_stat = -2*(-null_lklh + alt_lklh)
    
    #pvalue = 1 - scipy.stats.chi2.cdf(test_stat, 3)
    #print(pvalue)
    return test_stat


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
    all_teststats = []
    for i in range(num_snps):
        all_values = pvalues[i][1:]
        all_values = -np.log(all_values)
        meta_teststat = likelihood_ratio_test(all_values)
        all_teststats.append(meta_teststat)
    snps = pd.read_csv(input, usecols=[0], sep='\t', names = ['snp'], skiprows = 1)
    snps.insert(column = 'mixture_test_stats', loc = 1, value = all_teststats)
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
