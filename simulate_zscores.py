import argparse
import numpy as np
from numpy import linalg as LA
from numpy import genfromtxt

def simulateZscores(zfile, efile, qfile, output, n):
    #mean_zscores = np.loadtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_meanzscores.csv', dtype=complex, delimiter='\t')
    mean_zscores = np.loadtxt(zfile, dtype=complex, delimiter='\t')
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
    sim_zscores = []
    e_matrix = np.dot(Q, E)
    for i in range(n):
        if (i%10==0):
            print(i)
        z = np.random.normal(0, 1, n_genes)
        cur_sim_zscores = mean_zscores + np.dot(e_matrix,z)
        sim_zscores.append(cur_sim_zscores.real)
    print('simuated zscores calculated')

    sim_zscores = np.array(sim_zscores)
   # np.savetxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_simzscores100.csv', sim_zscores, delimiter='\t')
    np.savetxt(output, sim_zscores, delimiter='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-z", "--mzscores", required=True, help="Input mean zscores file")
    parser.add_argument("-e", "--eigenvalues", required=True, help="Input eigenvalues file")
    parser.add_argument("-q", "--eigenvectors", required=True, help="Input eigenvectorsfile")
    parser.add_argument("-o", "--output", required=True, help="Ouptput file with simulated zscores")
    parser.add_argument("-n", "--simulations", required=True, type=int, help="Number of simulations")
    params = parser.parse_args()
    simulateZscores(params.mzscores, params.eigenvalues, params.eigenvectors, params.output, params.simulations)


if __name__ == "__main__":
    main()
