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
from scipy import stats, linalg as LA
import scipy.optimize
from scipy.stats import chi2, norm
import sys
import multiprocessing as mp
from xQTL import __version__

def ERROR(msg):
    sys.stderr.write("[ERROR]: " + msg.strip() + "\n")
    sys.exit(1)

def MSG(msg):
    sys.stderr.write("[PROGRESS]: " + msg.strip() + "\n")

def eigendecomposition(sumstats):
    # read entire file into memory and obtain snp-gene matrix of t-stats
    MSG("Starting eigendecomposition null method")
    data = pd.read_csv(sumstats, sep='\t')
    zscores = data.pivot_table(index='SNP', columns='gene', values='t-stat')
    data = ''
    # get cov
    genes = list(zscores.columns)
    scores = np.transpose(np.array(zscores))
    zscores = ''
    cov = np.cov(scores)
    # get mean zscores
    mean_zscores = []
    for i in scores:
        mean_zscores.append(np.mean(i))
    mean_zscores = np.array(mean_zscores)
    # perform eigendecomposition
    e_values, Q = LA.eigh(cov)
    e_values = e_values.real
    e_values[e_values < 0] = 0
    Q = Q.real
    # get e_matrix
    n_genes = len(e_values)
    diag_e_values = np.diag(e_values)
    E = np.sqrt(diag_e_values)
    e_matrix = np.dot(Q, E)
    MSG("Finished eigendecomposition")
    return mean_zscores, e_matrix

def simulate_null_values(mean_zscores, e_matrix, cpma, xqtl, job_queue, writer_queue):
    MSG("Starting simulation of null values")
    n_genes = len(mean_zscores)
    # simulate 500000 null values
    num_sim = 500000
    # perform in chunks of 30000 for faster operations
    # can increase chunk size if more memory available
    iterations = math.ceil(num_sim/30000)
    sim_undone = num_sim
    index = 0
    for i in range(iterations):
        cur_n = min(30000, sim_undone)
        sim_undone = sim_undone - cur_n
        MSG(f"Null value simulation #{num_sim-sim_undone}")
        z=np.random.normal(0, 1, (n_genes, cur_n))
        mzscores_tile = np.transpose(np.tile(mean_zscores, (cur_n, 1)))
        sim_zscores = mzscores_tile + np.dot(e_matrix, z)
        for sim in np.transpose(sim_zscores):
            job_queue.put([f'snp_{index}', sim, cpma, xqtl])
            index += 1

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

    return test_stat, bestT, bestL

def GetPvalue(null_method, null_values, test_stat):
    if null_method == "chi2":
        pval = 1 - chi2.cdf(test_stat, 1)
    else:
        n_iter = len(null_values)
        index = n_iter - np.searchsorted(null_values, test_stat)
        pval = (index + 1)/(n_iter+1)
    return pval

def RunSNP(snp, tstats, CPMA=False, XQTL=False, \
        null_method="chi2", null_values_cpma=None, null_values_xqtl=None, null_sim=False):
    results = {}
    results["snp"] = snp
    pvalues = 2*stats.norm.cdf(-np.abs(tstats))
    if CPMA:
        cpma = calculate_cpma(pvalues)
        results["CPMA"] = cpma
        if not null_sim:
            results["CPMA_p"] = GetPvalue(null_method, null_values_cpma, cpma)
    if XQTL:
        test_stat, best_T, best_L = calculate_xqtl(pvalues)
        results["xQTL"] = test_stat
        results["predicted_T"] = best_T
        results["predicted_L"] = best_L
        if not null_sim:
            results["xQTL_p"] = GetPvalue(null_method, null_values_xqtl, test_stat)
    return results

def WriteSNP(results, outf, CPMA=False, XQTL=False, null_sim=False):
    outitems = [results["snp"]]
    if CPMA:
        outitems.extend([results["CPMA"]])
        if not null_sim:
            outitems.extend([results["CPMA_p"]])
    if XQTL:
        outitems.extend([results["xQTL"],results["predicted_T"], \
            results["predicted_L"]])
        if not null_sim:
            outitems.extend([results["xQTL_p"]])
    outf.write("\t".join([str(item) for item in outitems])+"\n")

def worker(job_queue, out_queue, null_method, null_values_cpma, null_values_xqtl, null_sim):
    # Keep looking for jobs to process. add results to out_queue
    while True:
        item = job_queue.get() # [snp, tstats, args.cpma, args.xqtl]
        if item == "DONE": break
        results = RunSNP(*item[0:4], null_method, \
            null_values_cpma, null_values_xqtl, null_sim)
        out_queue.put([results, item[2], item[3]])

def writer(out_queue, CPMA, XQTL, out, null_sim):
    # Set up output file
    header = ["SNP"]
    if CPMA:
        header.extend(["CPMA"])
        if not null_sim:
            header.extend(["CPMA_p"])
    if XQTL:
        header.extend(["xQTL","predicted_T","predicted_L"])
        if not null_sim:
            header.extend(["xQTL_p"])
    if not null_sim:
        outf = open(f'{out}.tab', "w")
    else:
        outf = open(f'{out}_null_sim.tab', "w")
    outf.write("\t".join(header)+"\n")

    # Keep looking for output jobs from the queue
    while True:
        item = out_queue.get() # [results, CPMA, XQTL]
        if item == "DONE": break
        WriteSNP(item[0], outf, item[1], item[2], null_sim)
    outf.close()

NULLOPTIONS = ["chi2","eigen"]
def main(args):
    if not os.path.exists(args.input):
        ERROR("Could not find %s"%args.input)
    if args.null_method not in NULLOPTIONS:
        ERROR("--null_method must be one of %s"%NULLOPTIONS)

    # if eigen method to get empirical null, set up job queue and processors to use in simulating null values
    if args.null_method == 'eigen':
        null_sim = True
        job_queue = mp.Queue()
        out_queue = mp.Queue()
        processes = [mp.Process(target=worker, args=(job_queue, out_queue, \
            args.null_method, None, None, null_sim)) \
            for i in range(np.max([args.threads-1, 1]))]
        writer_proc = mp.Process(target=writer, args=(out_queue, args.cpma, args.xqtl, args.out, null_sim))
        for p in processes: p.start()
        writer_proc.start()
        
        # get eigendecomposition and use job queue to simulate null values
        mean_zscores, e_matrix = eigendecomposition(args.input)
        simulate_null_values(mean_zscores, e_matrix, args.cpma, args.xqtl, job_queue, out_queue)
        
        for i in range(args.threads): job_queue.put("DONE")
        for p in processes: p.join()
        out_queue.put("DONE")
        writer_proc.join() 
        MSG("Finished getting empirical null")

    null_values_cpma = []
    null_values_xqtl = [] 
    with open(f'{args.out}_null_sim.tab') as f:
        header = f.readline()
        header_items = header.split()
        if args.cpma:
            cpma_col = header_items.index("CPMA")
        if args.xqtl:
            xqtl_col = header_items.index("xQTL")
        for line in f:
            values = line.strip().split('\t')
            if args.cpma:
                null_values_cpma.append(float(values[cpma_col]))
            if args.xqtl:
                null_values_xqtl.append(float(values[xqtl_col]))
    null_values_cpma.sort()
    null_values_xqtl.sort()

    # Set up job queue and processors
    # Note right now will use at least two processors
    # should possibly not do all this if threads==1
    null_sim = False
    job_queue = mp.Queue()
    out_queue = mp.Queue()
    processes = [mp.Process(target=worker, args=(job_queue, out_queue, \
        args.null_method, null_values_cpma, null_values_xqtl, null_sim)) \
        for i in range(np.max([args.threads-1, 1]))]
    writer_proc = mp.Process(target=writer, args=(out_queue, args.cpma, args.xqtl, args.out, null_sim))
    for p in processes: p.start()
    writer_proc.start()

    MSG("Starting to read user input")
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
    MSG("Finished")
    return 0

def getargs(): 
    parser = argparse.ArgumentParser(__doc__)
    run_group = parser.add_argument_group("xQTL run parameters")
    run_group.add_argument("--xqtl", help="Run x-QTL", action="store_true")
    run_group.add_argument("--cpma", help="Run CPMA", action="store_true")
#    run_group.add_argument("--eigendecomp", help="Run eigendecomposition method to adjust for gene correlation", action="store_true")
#    run_group.add_argument("--snpermute", help="Run snp-permute method to adjust for gene correlation", action="store_true")
#    run_group.add_argument("-s", "--seed", help="Seed for random generator", type=int, default=0 )    
    run_group.add_argument("--null_method", help="How to get the null distribution for test stats" \
            "Options: %s"%NULLOPTIONS, type=str, default="chi2")
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
