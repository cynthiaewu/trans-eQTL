import pandas as pd
import numpy as np
import scipy.optimize, scipy.stats
from numpy import genfromtxt
import math
import argparse


def calculate_cpma_topx(pvalues, x):
    #pvalues = genfromtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_pvalues_nofdr.csv', delimiter='\t', skip_header=1)
    #pvalues = genfromtxt(input, delimiter='\t', skip_header=1)
    #np.negative(np.log(pvalues[1:]))
    if len(pvalues.shape) == 1:
      pvalues = pvalues.reshape(1, -1)
    #print(len(pvalues.shape))
    num_snps = len(pvalues)
    #num_genes = len(pvalues[1])-1
    num_genes = len(pvalues[0])-1
    top_numgenes = math.ceil(x*num_genes)
    pvalues.sort()
    print(pvalues[0])
    all_values = -np.log(pvalues[0])
    top_values = all_values[:top_numgenes]
    bordernum = top_values[-1]
    #print(top_values)
    def log_likelihood_neg_cpma(x):
        #print(x)
        log_likelihood = top_numgenes*np.log(x) - \
                         x*np.sum(top_values) + \
                         (num_genes-top_numgenes)*np.log(1-np.exp(-x*bordernum))
        return -log_likelihood
    if math.isclose(x, 1.0):
        max_likelihood = 1/np.mean(all_values)
    else:
        max_likelihood = scipy.optimize.minimize(log_likelihood_neg_cpma,
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
    #return cpmax_value
    return pvalue
    #print(cpma)
    #snps = pd.read_csv(input, usecols=[0], sep='\t', names = ['snp'], skiprows = 1)
    #snps.insert(column = 'cpma', loc = 1, value = cpma)
    #snps.to_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_cpma_nofdr.csv', index=False, header=True, sep='\t')
    #snps.to_csv(output, index=False, header=True, sep='\t')



def log_likelihood_neg(t, L, pvals):
    pvals = np.array(pvals)
    return -np.sum(np.log((1 - t) * np.exp(-pvals) + t * 1/L * np.exp(-1/L*pvals)))


def likelihood_ratio_test(pvals, trueT=None, trueL=None):
    assert (trueT is None) == (trueL is None)
    null_lklh = log_likelihood_neg(0, 1, pvals)
    if trueT is None:
        alt_lklh = scipy.optimize.minimize(lambda tL: log_likelihood_neg(*tL, pvals=pvals),
                                           method='L-BFGS-B',
                                           x0=(0.5, 1),
                                           bounds=((10**(-5), 1),(10**(-5), None)))['fun']
    else:
        alt_lklh = log_likelihood_neg(t=trueT, L=trueL, pvals=pvals)
    test_stat = -2*(-null_lklh + alt_lklh)
    
    pvalue = 1 - scipy.stats.chi2.cdf(test_stat, 1)
    return pvalue


def get_bestTL_guesses(pvals, init_guesses=None):
    if init_guesses is None:
        init_guesses = np.linspace(0, 1, 101)[1:-1]
    max_lklh = float('inf')
    for guess in init_guesses:
        # print(guess)
        results = scipy.optimize.minimize(lambda tL: log_likelihood_neg(*tL, pvals=pvals),
                                           # method='L-BFGS-B',
                                           x0=(guess, 1),
                                           bounds=((10**(-5), 1-10**(-5)),(10**(-5), None)))
        if results['fun'] < max_lklh:
            max_lklh = results['fun']
            x = results['x']
            bestT, bestL = x[0], x[1]
        # print(f'{max_lklh} bestT, bestL: {bestT}, {bestL}')
    return bestT, bestL


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
    model = []
    inferredT = []
    for i in range(num_snps):
        all_values = pvalues[i][1:]
        pvals = -np.log(all_values)
        bestT, bestL = get_bestTL_guesses(pvals)
        if bestT > 0.05:
            #meta_pval = likelihood_ratio_test(all_values, bestT, bestL)
            meta_pval = likelihood_ratio_test(all_values)
            model.append('mixture')
            inferredT.append(-1)
            print(f'bestT: {bestT}, mixture model chosen')
        else:
            meta_pval = calculate_cpma_topx(all_values, bestT)
            model.append('cpmax')
            inferredT.append(bestT)
            print(f'bestT: {bestT}, cpmax chosen')
        all_metapvals.append(meta_pval)
    snps = pd.read_csv(input, usecols=[0], sep='\t', names = ['snp'], skiprows = 1)
    snps.insert(column = 'mixture_cpmax_pval', loc = 1, value = all_metapvals)
    snps.insert(column = 'model', loc = 2, value = model)
    snps.insert(column = 'inferredT', loc = 3, value = inferredT)
    #snps.to_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_cpma_nofdr.csv', index=False, header=True, sep='\t')
    #snps.to_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_cpma_nofdr.csv', index=False, header=True, sep='\t')
    snps.to_csv(output, index=False, header=True, sep='\t')
    


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input pvalues no fdr file")
    parser.add_argument("-o", "--output", required=True, help="Ouptput file with pvalues from testing if from mixtureModel if inferred T is > 0.05, else from cpmax with x as inferred T")
    params = parser.parse_args()
    get_pvalue(params.input, params.output)


if __name__ == "__main__":
    main()
