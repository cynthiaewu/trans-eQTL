import pandas as pd
import numpy as np
import math
from scipy.stats import norm
import argparse


def calculate_cpma(sim_pvalues, num_genes):
    likelihood = 1/(np.mean(np.negative(np.log(sim_pvalues))))
    value = -2 * ((((likelihood - 1) * num_genes)/likelihood) - num_genes*np.log(likelihood))
    return value


def sim_cpma_values(n, num_genes):
    iterations = math.ceil(n/50000)
    sim_undone = n
    sim_cpma = []
    for i in range(iterations):
        cur_n = min(50000, sim_undone)
        print(n-sim_undone)
        sim_undone = sim_undone - cur_n
        sim_zscores = np.random.normal(0, 1, (cur_n, num_genes))
        sim_pvalues = norm.cdf(sim_zscores)
        #print(sim_pvalues.shape)
        for sim in sim_pvalues:
            cpma = calculate_cpma(sim, num_genes)
            sim_cpma.append(cpma)
    return sim_cpma
 

def calculate_empirical(num_simvalues, num_genes, obs_file, output):
    #sim_zscores = np.loadtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_simzscores100.csv', delimiter='\t')
    cpma = sim_cpma_values(num_simvalues, num_genes)
    cpma.sort()
    #obs_cpma = pd.read_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_cpma_nofdr.csv', sep='\t')
    obs_cpma = pd.read_csv(obs_file, sep='\t')
    n_iter = len(cpma)
    pvalue = []
    for obs in obs_cpma['cpma']:
        index = n_iter - np.searchsorted(cpma, obs)
        pvalue.append((index + 1)/(n_iter+1))
    obs_cpma.insert(2, "pvalue", pvalue)
    
    obs_cpma.to_csv(output, index=False, header=True, sep='\t')
    print(f'Finished calculating empirical pvalues for {obs_file}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--num_simvalues", required=True, type=int, default=500000, help="Number of simulated cpma values for empirical null")
    parser.add_argument("-g", "--num_genes", required=True, type=int, default=15000, help="Number of genes")
    #parser.add_argument("-c", "--config", required=True, help="Input config file with parameters")
    parser.add_argument("-o", "--observed", required=True, help="Input observed cpma for each snp file")
    parser.add_argument("-e", "--emp_output", required=True, help="Ouptput file with empirical pvalues for each snp")
    params = parser.parse_args()
    calculate_empirical(params.num_simvalues, params.num_genes,  params.observed, params.emp_output)


if __name__ == "__main__":
    main()
