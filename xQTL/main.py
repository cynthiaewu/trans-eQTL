#!/usr/bin/env python3

"""
Tool to simulate expression data

Example command:
xQTL-run --input ../trans-eQTL/test/sim-0/matrixEQTL.out.gz --out test --cpma --xqtl 
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
import multiprocessing as mp
from xQTL import __version__

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

def RunSNP(snp, tstats, CPMA=False, XQTL=False):
    results = {}
    results["snp"] = snp
    pvalues = 2*stats.norm.cdf(-np.abs(tstats))
    if CPMA:
        cpma = calculate_cpma(pvalues)
        results["CPMA"] = cpma
        results["CPMA_p"] = None # TODO
    if XQTL:
        test_stat, best_T, best_L = calculate_xqtl(pvalues)
        results["xQTL"] = test_stat
        results["predicted_T"] = best_T
        results["predicted_L"] = best_L
        results["xQTL_p"] = None # TODO
    return results

def WriteSNP(results, outf, CPMA=False, XQTL=False):
    outitems = [results["snp"]]
    if CPMA:
        outitems.extend([results["CPMA"],results["CPMA_p"]])
    if XQTL:
        outitems.extend([results["xQTL"],results["predicted_T"], \
            results["predicted_L"], results["xQTL_p"]])
    outf.write("\t".join([str(item) for item in outitems])+"\n")

def worker(job_queue, out_queue):
    # Keep looking for jobs to process. add results to out_queue
    while True:
        item = job_queue.get() # [snp, tstats, args.cpma, args.xqtl]
        if item == "DONE": break
        results = RunSNP(*item[0:4])
        out_queue.put([results, item[2], item[3]])

def writer(out_queue, CPMA, XQTL, out):
    # Set up output file
    header = ["SNP"]
    if CPMA:
        header.extend(["CPMA","CPMA_p"])
    if XQTL:
        header.extend(["xQTL","predicted_T","predicted_L","xQTL_p"])
    outf = open(f'{out}.tab', "w")
    outf.write("\t".join(header)+"\n")

    # Keep looking for output jobs from the queue
    while True:
        item = out_queue.get() # [results, CPMA, XQTL]
        if item == "DONE": break
        WriteSNP(item[0], outf, item[1], item[2])

    outf.close()

def main(args):
    if not os.path.exists(args.input):
        ERROR("Could not find %s"%args.input)

    # Set up job queue and processors
    # Note right now will use at least two processors
    # should possibly not do all this if threads==1
    job_queue = mp.Queue()
    out_queue = mp.Queue()
    processes = [mp.Process(target=worker, args=(job_queue, out_queue)) \
        for i in range(np.max([args.threads-1, 1]))]
    writer_proc = mp.Process(target=writer, args=(out_queue, args.cpma, args.xqtl, args.out))
    for p in processes: p.start()
    writer_proc.start()

    MSG("Starting to read input")
    # TODO add check for SNP sorted
    with gzip.open(args.input,'rt') as f:
        header = f.readline() # SNP gene    beta    t-stat  p-value FDR
        header_items = header.split()
        snp_col = header_items.index("SNP")
        tstat_col = header_items.index("t-stat")
        snp = ''
        tstats = []
        for line in f:
            values = line.strip().split('\t')
            cur_snp = values[snp_col]
            if snp == cur_snp: 
                tstats.append(float(values[tstat_col]))
            elif snp:
                job_queue.put([snp, tstats, args.cpma, args.xqtl])
            snp = cur_snp
        # Make sure final SNP gets run
        if snp:
            job_queue.put([snp, tstats, args.cpma, args.xqtl])
    for i in range(args.threads): job_queue.put("DONE")
    for p in processes: p.join()
    out_queue.put("DONE")
    writer_proc.join()
    return 0

def getargs(): 
    parser = argparse.ArgumentParser(__doc__)
    run_group = parser.add_argument_group("xQTL run parameters")
    run_group.add_argument("--xqtl", help="Run x-QTL", action="store_true")
    run_group.add_argument("--cpma", help="Run CPMA", action="store_true")
#    run_group.add_argument("--eigendecomp", help="Run eigendecomposition method to adjust for gene correlation", action="store_true")
#    run_group.add_argument("--snpermute", help="Run snp-permute method to adjust for gene correlation", action="store_true")
#    run_group.add_argument("-s", "--seed", help="Seed for random generator", type=int, default=0 )    
    run_group.add_argument("--threads", help="Number of threads to use", type=int, default=1)
    inout_group = parser.add_argument_group("Input/output")
    inout_group.add_argument("--input", help="Matrix eQTL file. Must be sorted by SNP.", type=str, required=True)
    inout_group.add_argument("--out", help="Output prefix", type=str, required=True)
    ver_group = parser.add_argument_group("Version")
    ver_group.add_argument("--version", action="version", version = '{version}'.format(version=__version__))
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
