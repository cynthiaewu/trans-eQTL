import pandas as pd
import numpy as np
import argparse

def calculate_empirical(sim_file, obs_file, output, null):
    #sim_zscores = np.loadtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_simzscores100.csv', delimiter='\t')
    if not null:
        cpma = np.loadtxt(sim_file, delimiter='\t')
    else:
        cpma = np.loadtxt(sim_file, skiprows=1, usecols=1)
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
    parser.add_argument("-s", "--simulated", required=True, help="Input simulated cpma file")
    parser.add_argument("-o", "--observed", required=True, help="Input observed cpma for each snp file")
    parser.add_argument("-e", "--emp_output", required=True, help="Ouptput file with empirical pvalues for each snp")
    parser.add_argument("-n", "--null", default=False, help="If input simulated cpma file is from permuted null distribution")
    params = parser.parse_args()
    calculate_empirical(params.simulated, params.observed, params.emp_output, params.null)


if __name__ == "__main__":
    main()
