import pandas as pd
import numpy as np
import argparse

def calculate_empirical(null_teststats_file, obs_file, output):
    #sim_zscores = np.loadtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_simzscores100.csv', delimiter='\t')
    null_teststats = np.loadtxt(null_teststats_file, delimiter='\t')
    null_teststats.sort()
    #obs_cpma = pd.read_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_cpma_nofdr.csv', sep='\t')
    obs_teststats = pd.read_csv(obs_file, sep='\t')
    n_iter = len(null_teststats)
    pvalue = []
    for obs in obs_teststats['mixture_test_stats']:
        index = n_iter - np.searchsorted(null_teststats, obs)
        pvalue.append((index + 1)/(n_iter+1))
    obs_teststats.insert(2, "empirical_pvalue", pvalue)
    
    obs_teststats.to_csv(output, index=False, header=True, sep='\t')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--null", required=True, help="Input simulated null test stats file")
    parser.add_argument("-o", "--observed", required=True, help="Input observed test stats for each snp file")
    parser.add_argument("-e", "--emp_output", required=True, help="Ouptput file with empirical pvalues for each snp")
    params = parser.parse_args()
    calculate_empirical(params.null, params.observed, params.emp_output)


if __name__ == "__main__":
    main()
