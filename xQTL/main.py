#!/usr/bin/env python3

"""
Tool to simulate expression data

Example command:
xQTL-run --input test/sim-0/matrixEQTL.out --out test
"""

import argparse
import numpy as np
import os
import gzip
import pandas as pd
import math
from scipy import stats
import scipy.optimize
import sys

#from xQTL import __version__

def ERROR(msg):
    sys.stderr.write("[ERROR]: " + msg.strip() + "\n")
    sys.exit(1)

def MSG(msg):
    sys.stderr.write("[PROGRESS]: " + msg.strip() + "\n")

def calculate_cpma(pvalues):
    num_genes = len(pvalues)
    likelihood = 1/(np.mean(np.negative(np.log(pvalues))))
    value = -2 * ((((likelihood - 1) * num_genes)/likelihood) - num_genes*np.log(likelihood))
    if math.isnan(value):
        value = 0
    return value

def calculate_xqtl(pvalues):
    pvalues = np.where(pvalues > 10**(-150), pvalues, 10**(-150))
    values = -np.log(pvalues)
    test_stat, best_T, best_L = likelihood_ratio_test(values)
    return test_stat, best_T, best_L

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

# calculate pvalues based off of chi square distribution if there is no gene correlation
def calculate_pvalue_chidist_nocorr():
    return values

# simulate null distribution of xQTL or CPMA values to adjust for gene correlation
# permute genotypes of existing snps randomly to create null snps and obtain xQTL or CPMA null values
def simulate_snp_permute():
    return values

# simulate null distribution of xQTL or CPMA values to adjust for gene correlation
# perform eigendecomposition to use decomposed elements to simulate xQTL or CPMA null values
def simulate_eigendecomp():
    return values

# calculate empirical pvalues with simulated xQTL or CPMA null distribution
def calculate_empirical_pvalue():
    return values

def main(args):
    if not os.path.exists(args.input):
        ERROR("Could not find %s"%args.input)

    # Load Matrix eQTL data
    # TODO - require matrix eQTL results sorted by SNP, then analyze data for one SNP at a time
    # can parralelize by SNP
    #data = pd.read_csv(args.input, sep="\t")
    if args.cpma:
        cpma_out = open(f'{args.out}_cpma', "w")
        cpma_out.write(f'SNP\tCPMA\n')
    if args.xqtl:
        xqtl_out = open(f'{args.out}_xqtl', "w")
        xqtl_out.write(f'SNP\txQTL\tpredicted_T\tpredicted_L\n')

    MSG("Starting to read input")
    with gzip.open(args.input,'rt') as f:
        #header = f.readline().decode("utf-8")
        header = f.readline()
        snp = ''
        tstats = []
        for line in f:
            #line = line.decode("utf-8")
            values = line.strip().split('\t')
            cur_snp = values[0]
            if snp == cur_snp: 
                tstats.append(float(values[3]))
            elif snp:
                pvalues = 2*stats.norm.cdf(-np.abs(tstats))
                if args.cpma:
                    cpma = calculate_cpma(pvalues)
                    cpma_out.write(f'{snp}\t{cpma}\n')
                    # write to output
                if args.xqtl:
                    test_stat, best_T, best_L = calculate_xqtl(pvalues)
                    xqtl_out.write(f'{snp}\t{test_stat}\t{best_T}\t{best_L}\n')
            snp = cur_snp
    if args.cpma:
        cpma_out.close()
    if args.xqtl:
        xqtl_out.close()
                
    #zscores = data.pivot_table(index='SNP', columns='gene', values='t-stat')

    # Convert tvals to pvals
    #print(zscores) # TODO figure this out
    #pvalues = 2*stats.norm.cdf(-np.abs(zscores))
    #print(pvalues)
    #zscores["pval"] = 2*stats.norm.cdf(-np.abs(zscores))

    # TODO - compute test statistics, either 
    # using CPMA  ./calculate_cpma.py
    # or xQTL ./calculate_mixturemodel.py
    # Add options for either

    # TODO - get pval on test stat

    #ERROR("Not implemented")

def getargs(): 
    parser = argparse.ArgumentParser(__doc__)
    run_group = parser.add_argument_group("xQTL run parameters")
    run_group.add_argument("--xqtl", help="Run x-QTL", action="store_true")
    run_group.add_argument("--cpma", help="Run CPMA", action="store_true")
    run_group.add_argument("--eigendecomp", help="Run eigendecomposition method to adjust for gene correlation", action="store_true")
    run_group.add_argument("--snpermute", help="Run snp-permute method to adjust for gene correlation", action="store_true")
    run_group.add_argument("-s", "--seed", help="Seed for random generator", type=int, default=0 )    
    inout_group = parser.add_argument_group("Input/output")
    inout_group.add_argument("--input", help="Matrix eQTL file", type=str, required=True)
    inout_group.add_argument("--out", help="Output prefix", type=str, required=True)
    #ver_group = parser.add_argument_group("Version")
    #ver_group.add_argument("--version", action="version", version = '{version}'.format(version=__version__))
    args = parser.parse_args()
    return args

def run():
    args = getargs()
    if args == None:
        sys.exit(1)
    else:
        retcode = main(args)
        sys.exit(retcode)

if __name__ == "__main__":
    run()
