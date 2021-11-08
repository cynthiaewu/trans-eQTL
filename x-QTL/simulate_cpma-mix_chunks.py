import argparse
import numpy as np
from numpy import linalg as LA
from scipy.stats import norm
import scipy.optimize, scipy.stats
import math


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
    return test_stat


def calculate_cpma(sim_pvalues, num_genes):
    likelihood = 1/(np.mean(np.negative(np.log(sim_pvalues))))
    value = -2 * ((((likelihood - 1) * num_genes)/likelihood) - num_genes*np.log(likelihood))
    return value


def simulateZscores(zfile, efile, qfile, output_cpma, output_cpmamix, n):
    #mean_zscores = np.loadtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_meanzscores.csv', dtype=complex, delimiter='\t')
    mean_zscores = np.loadtxt(zfile, delimiter='\t')
    print(f'mean zscores file read {zfile}')
    #print(mean_zscores.shape)

    #e_values = np.loadtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_evalues.csv', dtype=complex, delimiter='\t')
    #e_values = (np.loadtxt(efile, dtype=complex, delimiter='\t')).real
    e_values = (np.loadtxt(efile, delimiter='\t'))
    n_genes = len(e_values)
    print(n_genes)
    print(f'e_values file read {efile}')
    #Q = np.loadtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_Q.csv', dtype=complex, delimiter='\t')
    #Q = (np.loadtxt(qfile, dtype=complex, delimiter='\t')).real
    Q = (np.loadtxt(qfile, delimiter='\t'))
    print(f'Q file read {qfile}')
    diag_e_values = np.diag(e_values)
    E = np.sqrt(diag_e_values)
    #print(E)
    
    print(f'starting simulations {zfile}')
    #print(Q.shape)
    #print(E.shape)
    e_matrix = np.dot(Q, E)
    Q = ''
    E = ''
   
    sim_cpma = []
    sim_cpmamix = []
    #iterations = math.ceil(n/30000)
    iterations = math.ceil(n/30000)
    #print(iterations)
    sim_undone = n
    #perform in chunks of 1000
    for i in range(iterations):
        #cur_n = min(30000, sim_undone)
        cur_n = min(30000, sim_undone)
        sim_undone = sim_undone - cur_n
        print(f'{n-sim_undone} {zfile}')

        z=np.random.normal(0, 1, (n_genes, cur_n))
        #print(z.shape)
        mzscores_tile = np.transpose(np.tile(mean_zscores, (cur_n, 1)))
        #print(mzscores_tile.shape)
        sim_zscores = mzscores_tile + np.dot(e_matrix, z)
        #print(sim_zscores.shape)
        sim_pvalues = np.transpose(2*norm.cdf(-np.abs(sim_zscores)))
        #sim_pvalues = np.transpose(norm.cdf(sim_zscores))
        #print(sim_pvalues.shape)
        for sim in sim_pvalues:
            #print(len(sim))
            cpmamix = likelihood_ratio_test(-np.log(sim))
            cpma = calculate_cpma(sim, n_genes)
            sim_cpma.append(cpma)
            sim_cpmamix.append(cpmamix)
        #print(f'iterations: {iterations}')

    print(f'simulated cpma calculated {zfile}')

    sim_cpma = np.array(sim_cpma)
    sim_cpmamix = np.array(sim_cpmamix)
   # np.savetxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_simzscores100.csv', sim_zscores, delimiter='\t')
    np.savetxt(output_cpma, sim_cpma, delimiter='\t', fmt='%f')
    np.savetxt(output_cpmamix, sim_cpmamix, delimiter='\t', fmt='%f')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-z", "--mzscores", required=True, help="Input mean zscores file")
    parser.add_argument("-e", "--eigenvalues", required=True, help="Input eigenvalues file")
    parser.add_argument("-q", "--eigenvectors", required=True, help="Input eigenvectorsfile")
    parser.add_argument("-o", "--output_cpma", required=True, help="Ouptput file with simulated cpma values")
    parser.add_argument("-m", "--output_cpmamix", required=True, help="Ouptput file with simulated cpma-mix values")
    parser.add_argument("-n", "--simulations", required=True, type=int, help="Number of simulations")
    params = parser.parse_args()
    np.random.seed(0)
    simulateZscores(params.mzscores, params.eigenvalues, params.eigenvectors, params.output_cpma, params.output_cpmamix, params.simulations)


if __name__ == "__main__":
    main()
