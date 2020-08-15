import argparse
import numpy as np
from numpy import linalg as LA
from scipy.stats import norm
#from scipy.stats import multivariate_normal
from numpy.random import standard_normal
from scipy.linalg import cholesky
import math


def calculate_cpma(sim_pvalues, num_genes):
    likelihood = np.mean(np.negative(np.log(sim_pvalues)))
    value = -2 * ((((likelihood - 1) * num_genes)/likelihood) - num_genes*np.log(likelihood))
    return value


def simulateZscores(zfile, covfile, output, n):
    #mean_zscores = np.loadtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_meanzscores.csv', dtype=complex, delimiter='\t')
    mean_zscores = np.loadtxt(zfile, delimiter='\t')
    print('mean zscores file read')
    #print(mean_zscores.shape)

    cov = np.loadtxt(covfile, delimiter='\t', skiprows=1)
    print('covariance file read')

    num_genes = len(mean_zscores)
    sim_cpma = []
    l = cholesky(cov, check_finite=False, overwrite_a=True)
    for i in range(n):
        sim_zscores = mean_zscores + l.dot(standard_normal(len(mean_zscores)))
        #sim_zscores = multivariate_normal(mean_zscores, cov)
        #print(sim_zscores)
        sim_pvalues = norm.cdf(sim_zscores)
        #print(sim_pvalues)
        cpma = calculate_cpma(sim_pvalues, num_genes)
        #print(cpma)
        sim_cpma.append(cpma)
        if i % 5000 == 0:
            print(f'Simulating {i}th cpma value')

    print('simulated cpma calculated')

    sim_cpma = np.array(sim_cpma)
   # np.savetxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_simzscores100.csv', sim_zscores, delimiter='\t')
    np.savetxt(output, sim_cpma, delimiter='\t', fmt='%f')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-z", "--mzscores", required=True, help="Input mean zscores file")
    parser.add_argument("-c", "--cov", required=True, help="Input covariance matrix file")
    parser.add_argument("-o", "--output", required=True, help="Ouptput file with simulated cpma values")
    parser.add_argument("-n", "--simulations", required=True, type=int, help="Number of simulations")
    params = parser.parse_args()
    np.random.seed(0)
    simulateZscores(params.mzscores, params.cov, params.output, params.simulations)


if __name__ == "__main__":
    main()
