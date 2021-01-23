import argparse
import numpy as np
from numpy import linalg as LA
from scipy.stats import norm


def calculate_cpma(sim_zscores, num_genes):
    sim_pvalues = norm.cdf(sim_zscores)
    likelihood = np.mean(np.negative(np.log(sim_pvalues)))
    value = -2 * ((((likelihood - 1) * num_genes)/likelihood) - num_genes*np.log(likelihood))
    return value


def simulateZscores(zfile, efile, qfile, output, n):
    #mean_zscores = np.loadtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_meanzscores.csv', dtype=complex, delimiter='\t')
    mean_zscores = np.loadtxt(zfile, delimiter='\t')
    print('mean zscores file read')

    #e_values = np.loadtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_evalues.csv', dtype=complex, delimiter='\t')
    e_values = np.loadtxt(efile, dtype=complex, delimiter='\t')
    n_genes = len(e_values)
    print(n_genes)
    print('e_values file read')
    #Q = np.loadtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_Q.csv', dtype=complex, delimiter='\t')
    Q = np.loadtxt(qfile, dtype=complex, delimiter='\t')
    print('Q file read')
    diag_e_values = np.diag(e_values)
    E = np.sqrt(diag_e_values)

    print('starting simulations')
    sim_cpma = []
    e_matrix = np.dot(Q, E)
    for i in range(n):
        if (i%1000==0):
            print(i)
        z = np.random.normal(0, 1, n_genes)
        sim_zscores = mean_zscores + np.dot(e_matrix,z)
        cpma = calculate_cpma(sim_zscores.real, n_genes)
        sim_cpma.append(cpma)
    print('simuated cpma calculated')

    sim_cpma = np.array(sim_cpma)
   # np.savetxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_simzscores100.csv', sim_zscores, delimiter='\t')
    np.savetxt(output, sim_cpma, delimiter='\t', fmt='%f')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-z", "--mzscores", required=True, help="Input mean zscores file")
    parser.add_argument("-e", "--eigenvalues", required=True, help="Input eigenvalues file")
    parser.add_argument("-q", "--eigenvectors", required=True, help="Input eigenvectorsfile")
    parser.add_argument("-o", "--output", required=True, help="Ouptput file with simulated cpma values")
    parser.add_argument("-n", "--simulations", required=True, type=int, help="Number of simulations")
    params = parser.parse_args()
    simulateZscores(params.mzscores, params.eigenvalues, params.eigenvectors, params.output, params.simulations)


if __name__ == "__main__":
    main()
