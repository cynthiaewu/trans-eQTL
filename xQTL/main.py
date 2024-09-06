#!/usr/bin/env python3

"""
Tool to simulate expression data

Example command:
xQTL-run --input ../trans-eQTL/test/sim-0/matrixEQTL.out.gz --out test --cpma --xqtl 
"""

import argparse
import warnings
import numpy as np
import os
import gzip
import pandas as pd
import pathlib
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
    num_snps = len(zscores.index)
    if num_snps < 2:
        ERROR("Input file has less than 2 snps, not enough for eigendecomposition null method")
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

def simulate_null_values(mean_zscores, e_matrix, cpma, xqtl, grid, job_queue, writer_queue):
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
        z=np.random.normal(0, 1, (n_genes, cur_n))
        mzscores_tile = np.transpose(np.tile(mean_zscores, (cur_n, 1)))
        sim_zscores = mzscores_tile + np.dot(e_matrix, z)
        for sim in np.transpose(sim_zscores):
            job_queue.put([f'snp_{index}', sim, cpma, xqtl, grid])
            index += 1

def calculate_cpma(pvalues):
    pvalues = np.where(pvalues > 10**(-150), pvalues, 10**(-150))
    num_genes = len(pvalues)
    likelihood = 1/(np.mean(np.negative(np.log(pvalues))))
    value = -2 * ((((likelihood - 1) * num_genes)/likelihood) - num_genes*np.log(likelihood))
    return value

def calculate_xqtl(pvalues, grid):
    pvalues = np.where(pvalues > 10**(-150), pvalues, 10**(-150))
    values = -np.log(pvalues)
    test_stat, best_T, best_L = likelihood_ratio_test(values, grid)
    return test_stat, best_T, best_L

def log_likelihood_neg(t, L, pvals, C=1):
    pvals = np.array(pvals)
    return -np.sum(np.log((1 - t) * np.exp(-pvals) + t * 1/L * np.exp(-1/L*pvals)))-C*np.log(4*t*(1-t))
    #return -np.sum(np.log((1 - t) * np.exp(-pvals) + t * 1/L * np.exp(-1/L*pvals)))

def grid_scipy(pvals,
            nt=5, min_t=0.01, max_t=0.99,
            nL=5, min_L=0.1, max_L=100):
    ts = np.linspace(min_t, max_t, nt)
    Ls = np.linspace(min_L, max_L, nL)
    best_f = {'fun': np.inf}
    for t0 in ts:
        for L0 in Ls:
            res = scipy.optimize.minimize(lambda tL: log_likelihood_neg(*tL, pvals=pvals),
                        method='L-BFGS-B',
                        x0=(t0, L0),
                        bounds=((10**(-5), 1-10**(-5)),(10**(-5), None)))
            if (res['fun'] < best_f['fun']):
                best_f = res
    return best_f

def likelihood_ratio_test(pvals, grid, trueT=None, trueL=None):
    assert (trueT is None) == (trueL is None)
    null_lklh = log_likelihood_neg(0.5, 1, pvals)
    if trueT is None:
        if grid:
            results = grid_scipy(pvals)
        else:
            results = scipy.optimize.minimize(lambda tL: log_likelihood_neg(*tL, pvals=pvals),
                                           method='L-BFGS-B',
                                           x0=(0.5, 1),
                                           bounds=((10**(-5), 1-10**(-5)),(10**(-5), None)))
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

def GetDist(chrom, pos, gene, genes_dict):
    if chrom == genes_dict[gene][0]:
        return abs(pos-genes_dict[gene][1])
    else: return 10000000000

def RunSNP(snp, genes, tstats, genes_dict, n_genes, remove_closest=False, CPMA=False, XQTL=False, grid=False,\
        null_method="chi2", null_values_cpma=None, null_values_xqtl=None, null_sim=False):
    results = {}
    results["snp"] = snp
    #results["gc"] = gc
    if remove_closest: 
        chrom = snp.split('_')[0]
        pos = int(snp.split('_')[1])
        dist = [GetDist(chrom, pos, gene, genes_dict) for gene in genes]
        #n_min_values = 303
        n_min_values = n_genes
        indices = sorted(range(len(dist)), key=lambda k: dist[k])[:n_min_values]
        pvalues = np.array([element for i, element in enumerate(tstats) if i not in indices])
    else:
        pvalues = np.array(tstats)
    #tstats are genomic control adjusted pvalues instead
    #pvalues_tstats = 2*stats.norm.cdf(-np.abs(tstats))
    #pvalues = np.array([tstats[i] for i in indices])
    if CPMA:
        cpma = calculate_cpma(pvalues)
        results["CPMA"] = cpma
        if not null_sim:
            results["CPMA_p"] = GetPvalue(null_method, null_values_cpma, cpma)
    if XQTL:
        test_stat, best_T, best_L = calculate_xqtl(pvalues, grid)
        results["xQTL"] = test_stat
        results["predicted_T"] = best_T
        results["predicted_L"] = best_L
        if not null_sim:
            results["xQTL_p"] = GetPvalue(null_method, null_values_xqtl, test_stat)
    return results

def WriteSNP(results, outf, CPMA=False, XQTL=False, null_sim=False, null_values_cpma=None, null_values_xqtl=None):
    if not null_sim:
        outitems = [results["snp"]]
        #outitems.extend([results["gc"]])
    else:
        outitems = []
    if CPMA:
        outitems.extend([results["CPMA"]])
        if not null_sim:
            outitems.extend([results["CPMA_p"]])
        # save null values in memory in addition to writing it to file
        else:
            null_values_cpma.append([results["CPMA"]])
    if XQTL:
        outitems.extend([results["xQTL"]])
        if not null_sim:
            outitems.extend([results["predicted_T"], results["predicted_L"], results["xQTL_p"]])
        else:
            null_values_xqtl.append([results["xQTL"]])
    outf.write("\t".join([str(item) for item in outitems])+"\n")

def worker(job_queue, out_queue, null_method, null_values_cpma, null_values_xqtl, null_sim):
    # Keep looking for jobs to process. add results to out_queue
    while True:
        item = job_queue.get() # [snp, gc, genes, tstats, genes_dict, args.cpma, args.xqtl]
                #[snp, genes, tstats, genes_dict, args.n_genes, args.remove_closest, args.cpma, args.xqtl, args.grid])
        if item == "DONE": break
        #results = RunSNP(*item[0:8], null_method, \
        results = RunSNP(*item[0:9], null_method, \
            null_values_cpma, null_values_xqtl, null_sim)
        out_queue.put([results, item[6], item[7]])

def writer(out_queue, CPMA, XQTL, out, null_sim, null_values_cpma, null_values_xqtl):
    # Set up output file
    if not null_sim:
        header = ["SNP"]
        #header.extend(["gc"])
    else:
        header = []
    if CPMA:
        header.extend(["CPMA"])
        if not null_sim:
            header.extend(["CPMA_p"])
    if XQTL:
        header.extend(["xQTL"])
        if not null_sim:
            header.extend(["predicted_T","predicted_L","xQTL_p"])
    if not null_sim:
        outf = open(f'{out}/results.tab', "w")
    else:
        outf = open(f'{out}/null_sim.tab', "w")
    outf.write("\t".join(header)+"\n")
    
    counter = 0
    # Keep looking for output jobs from the queue
    while True:
        counter += 1
        if (counter % 50000 == 0) and null_sim:
            MSG(f"Null value simulation #{counter}")
        if (counter % 1000 == 0) and not null_sim:
            MSG(f"Finished {counter} SNPs")
        item = out_queue.get() # [results, CPMA, XQTL]
        if item == "DONE":
            if not null_sim:
                MSG(f"Finished {counter-1} SNPs")
            break
        WriteSNP(item[0], outf, item[1], item[2], null_sim, null_values_cpma, null_values_xqtl)
    outf.close()

NULLOPTIONS = ["chi2","eigen"]
def main(args):
    if not os.path.exists(args.input):
        ERROR("Could not find %s"%args.input)
    if args.null_method not in NULLOPTIONS:
        ERROR("--null_method must be one of %s"%NULLOPTIONS)
    
    # Set up output directory
    if not os.path.exists(args.out):
        os.mkdir(args.out)

    null_values_cpma = []
    null_values_xqtl = [] 
    # if eigen method to get empirical null, set up job queue and processors to use in simulating null values
    if args.null_method == 'eigen':
        # if user already has empirical null values file, read values from file
        # else simulate null values and write to file
        if args.precomputed_null:
        #if os.path.exists(f'{args.out}_null_sim.tab'):
            with open(f'{args.precomputed_null}') as f:
                MSG(f'Reading precomputed file for empirical null values')
                #MSG(f'Reading {args.out}_null_sim.tab for empirical null values')
                header = f.readline()
                header_items = header.split()
                if args.cpma:
                    #check if header items exist
                    if "CPMA" in header_items:
                        cpma_col = header_items.index("CPMA")
                    else:
                        ERROR(f"CPMA column does not exist in input {args.precomputed_null}")                   
                if args.xqtl:
                    if "xQTL" in header_items:
                        xqtl_col = header_items.index("xQTL")
                    else:
                        ERROR(f"xQTL column does not exist in input {args.precomputed_null}")                   
                for line in f:
                    values = line.strip().split('\t')
                    if args.cpma:
                        null_values_cpma.append(float(values[cpma_col]))
                    if args.xqtl:
                        null_values_xqtl.append(float(values[xqtl_col]))
        else:
            # get eigendecomposition
            np.random.seed(args.seed)
            mean_zscores, e_matrix = eigendecomposition(args.input)
            null_sim = True
            job_queue = mp.Queue()
            out_queue = mp.Queue()
            processes = [mp.Process(target=worker, args=(job_queue, out_queue, \
                args.null_method, None, None, null_sim)) \
                for i in range(np.max([args.threads-1, 1]))]
            writer_proc = mp.Process(target=writer, args=(out_queue, args.cpma, \
                args.xqtl, args.out, null_sim, null_values_cpma, null_values_xqtl))
            for p in processes: p.start()
            writer_proc.start()
            
            # use job queue to simulate null values
            simulate_null_values(mean_zscores, e_matrix, args.cpma, args.xqtl, args.grid, job_queue, out_queue)
            
            for i in range(args.threads): job_queue.put("DONE")
            for p in processes: p.join()
            out_queue.put("DONE")
            writer_proc.join() 
            MSG("Finished getting empirical null")

        null_values_cpma.sort()
        null_values_xqtl.sort()

    MSG("Starting to read user input")
    # TODO add check for SNP sorted
    input_fn = pathlib.Path(args.input)
    open_method = gzip.open if input_fn.suffix == '.gz' else open
    with open_method(input_fn,'rt') as f:
        header = f.readline() # SNP gene    beta    t-stat  p-value FDR
        header_items = header.split()
        #if not all(col in header_items for col in ["SNP", "gene", "t-stat"]): 
        if not all(col in header_items for col in ["SNP", "gene", "p-value"]): 
        #if not all(col in header_items for col in ["SNP", "gene", "p-value_gc"]): 
            ERROR(f"Required columns do not exist in input {args.input}")                   
        snp_col = header_items.index("SNP")
        gene_col = header_items.index("gene")
        #gc_col = header_items.index("gc_value")

        # use genomic control adjusted pvals instead of t-stats. t-stats were only used for eigendecomposition which we are not using
        #tstat_col = header_items.index("t-stat")
        #tstat_col = header_items.index("p-value_gc")
        tstat_col = header_items.index("p-value")
        # Set up job queue and processors
        # Note right now will use at least two processors
        # should possibly not do all this if threads==1
        null_sim = False
        job_queue = mp.Queue()
        out_queue = mp.Queue()
        processes = [mp.Process(target=worker, args=(job_queue, out_queue, \
            args.null_method, null_values_cpma, null_values_xqtl, null_sim)) \
            for i in range(np.max([args.threads-1, 1]))]
        writer_proc = mp.Process(target=writer, args=(out_queue, args.cpma, args.xqtl, args.out, null_sim, None, None))
        for p in processes: p.start()
        writer_proc.start()
        
        genes_dict = {}
        if args.remove_closest:
            genes_info = pd.read_csv(args.genes_info, sep="\t")
            for index, row in genes_info.iterrows():
                genes_dict[row['gene']] = [row['chrom'], row['avg.coord']]
    
        tstats = []
        genes = []
        # read the first snp and append the first tstat value
        firstsnp = f.readline()
        values = firstsnp.strip().split('\t')
        snp = values[snp_col]
        genes.append(values[gene_col])
        #gc = values[gc_col]
        tstats.append(float(values[tstat_col]))
        for line in f:
            values = line.strip().split('\t')
            cur_snp = values[snp_col]
            if snp == cur_snp: 
                genes.append(values[gene_col])
                tstats.append(float(values[tstat_col]))
            elif snp:
                #job_queue.put([snp, gc, genes, tstats, genes_dict, args.cpma, args.xqtl, args.grid])
                job_queue.put([snp, genes, tstats, genes_dict, args.n_genes, args.remove_closest, args.cpma, args.xqtl, args.grid])
                tstats = []
                genes = []
            snp = cur_snp
            #gc = values[gc_col]
        # Make sure final SNP gets run
        if snp:
            #job_queue.put([snp, gc, genes, tstats, genes_dict, args.cpma, args.xqtl, args.grid])
            job_queue.put([snp, genes, tstats, genes_dict, args.n_genes, args.remove_closest, args.cpma, args.xqtl, args.grid])
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
    run_group.add_argument("--grid", help="Run x-QTL with a grid search for t target genes", action="store_true")
    run_group.add_argument("--cpma", help="Run CPMA", action="store_true")
    run_group.add_argument("-s", "--seed", help="Seed for random generator", type=int, default=0 )    
    run_group.add_argument("--null_method", help="How to get the null distribution for test stats" \
            "Options: %s"%NULLOPTIONS, type=str, default="chi2")
    run_group.add_argument("--precomputed_null", help="Precomputed file with null values for test stats", \
            type=str)
    run_group.add_argument("--threads", help="Number of threads to use", type=int, default=1)
    run_group.add_argument("--remove_closest", help="Remove n closest genes for each variant", action="store_true")
    run_group.add_argument("--genes_info", help="File with info for genes. Must have 'chrom' and 'avg.coord'.", type=str)
    run_group.add_argument("--n_genes", help="Number of n closest genes to remove", type=int, default=1)
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
