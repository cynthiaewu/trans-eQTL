import numpy as np
from numpy import genfromtxt
import scipy.optimize, scipy.stats
import sys
import subprocess
import os
import argparse


def get_T(pvalue_file):
    #pvalues = genfromtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_pvalues_nofdr.csv', delimiter='\t', skip_header=1)
    pvalues = genfromtxt(pvalue_file, delimiter='\t', skip_header=1)
    #np.negative(np.log(pvalues[1:]))
    if len(pvalues.shape) == 1:
      pvalues = pvalues.reshape(1, -1)
    #print(len(pvalues.shape))
    num_snps = len(pvalues)
    #num_genes = len(pvalues[1])-1
    num_genes = len(pvalues[0])-1
    all_Ts = []
    for i in range(num_snps):
        all_values = pvalues[i][1:]
        all_values = -np.log(all_values)
        all_Ts.append(get_bestTL_guesses(all_values)[0])
    return all_Ts
        

def log_likelihood_neg(t, L, pvals):
    pvals = np.array(pvals)
    return -np.sum(np.log((1 - t) * np.exp(-pvals) + t * 1/L * np.exp(-1/L*pvals)))


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

'''
def get_bestTL(pvals):
    x = scipy.optimize.minimize(lambda tL: log_likelihood_neg(*tL, pvals=pvals),
                        method='L-BFGS-B',
                        x0=(0.5, 1),
                        bounds=((10**(-5), 1),(10**(-5), None)))['x']
    bestT, bestL = x[0], x[1]
    return bestT, bestL    
'''

def iterate_folders(folder, iterations):
    all_Ts_iter = []
    for i in range(iterations):
        input_folder = f'{folder}/Simulation_{i}'
        print(f'Starting mixture infer T for Simulation {i}, {folder}')
        cpma_folder = os.path.join(input_folder, 'CPMA')
        mixture_folder = os.path.join(input_folder, 'mixtureModel')
        if not os.path.isdir(mixture_folder):
            os.mkdir(mixture_folder)
        pvalue_file = f'{cpma_folder}/gene-snp-eqtl_pvalue'
        all_Ts = get_T(pvalue_file)
        all_Ts_iter.append(all_Ts)
    all_Ts_iter = np.array(all_Ts_iter).T
    np.savetxt(f'{folder}/inferred_mixtureT_guesses', all_Ts_iter)
   

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--folder', type=str, help='Input folder with simulated expression and genotype files')
    #parser.add_argument('-p', '--scripts_folder', type=str, help='Input folder with scripts')
    parser.add_argument('-i', '--iterations', type=int, help='# of simulated folders to run cpma pipeline on')
    params = parser.parse_args()

    iterate_folders(folder=params.folder, iterations=params.iterations)


if __name__ == '__main__':
    main()
