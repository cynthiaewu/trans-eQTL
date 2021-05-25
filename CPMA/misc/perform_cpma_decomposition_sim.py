import pandas as pd
import numpy as np
from scipy import linalg as LA
from scipy.stats import norm
import math
import argparse
import time


def getCov(input):
    data = pd.read_csv(input, sep='\t', index_col=0)
    genes = list(data.columns)
    scores = np.transpose(np.array(data))
    
    cov = np.cov(scores)
    #cov_matrix = pd.DataFrame(cov, columns=genes)
    #cov_matrix.to_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_cov.csv', index=True, header=True, sep='\t')
    #cov_matrix.to_csv(cov_out, index=False, header=True, sep='\t')
    print('cov matrix calculated')

    mean_zscores = []
    for i in scores:
        mean_zscores.append(np.mean(i))
    mean_zscores = np.array(mean_zscores)
    #np.savetxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_meanzscores.csv', mean_zscores, delimiter='\t')
    #np.savetxt(m_out, mean_zscores, delimiter='\t', fmt='%f')
    print('mean zscores calculated')
    return cov, mean_zscores


def calculate_values(cov):
    #zscores = genfromtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_zscores.csv', delimiter='\t', skip_header=1)
    #cov = np.loadtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_cov.csv', delimiter='\t', skiprows=1)
    #cov = np.loadtxt(covfile, delimiter='\t', skiprows=1)
    #print('files read')

    e_values, Q = LA.eig(cov)
    e_values = e_values.real
    e_values[e_values < 0] = 0
    Q = Q.real
    #np.savetxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_evalues.csv', e_values, delimiter='\t')
    #np.savetxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_Q.csv', Q, delimiter='\t')
    #np.savetxt(e_out, e_values, delimiter='\t')
    #np.savetxt(q_out, Q, delimiter='\t')
    #diag_e_values = np.diag(e_values)
    #E = np.sqrt(diag_e_values)
    print('eigendecomposition values calculated')
    return e_values, Q


def calculate_cpma(sim_pvalues, num_genes):
    likelihood = 1/(np.mean(np.negative(np.log(sim_pvalues))))
    value = -2 * ((((likelihood - 1) * num_genes)/likelihood) - num_genes*np.log(likelihood))
    return value


def simulate_cpma(mean_zscores, e_values, Q, output, n):
    #mean_zscores = np.loadtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_meanzscores.csv', dtype=complex, delimiter='\t')
    #mean_zscores = np.loadtxt(zfile, delimiter='\t')
    #print('mean zscores file read')
    #print(mean_zscores.shape)

    #e_values = np.loadtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_evalues.csv', dtype=complex, delimiter='\t')
    #e_values = (np.loadtxt(efile, dtype=complex, delimiter='\t')).real
    #e_values = (np.loadtxt(efile, delimiter='\t'))
    n_genes = len(e_values)
    print(n_genes)
    #print('e_values file read')
    #Q = np.loadtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_Q.csv', dtype=complex, delimiter='\t')
    #Q = (np.loadtxt(qfile, dtype=complex, delimiter='\t')).real
    #Q = (np.loadtxt(qfile, delimiter='\t'))
    #print('Q file read')
    diag_e_values = np.diag(e_values)
    E = np.sqrt(diag_e_values)
    #print(E)
    
    print('starting simulations')
    #print(Q.shape)
    #print(E.shape)
    e_matrix = np.dot(Q, E)
   
    sim_cpma = []
    #iterations = math.ceil(n/30000)
    iterations = math.ceil(n/30000)
    #print(iterations)
    sim_undone = n
    #perform in chunks of 1000
    for i in range(iterations):
        #cur_n = min(30000, sim_undone)
        cur_n = min(30000, sim_undone)
        sim_undone = sim_undone - cur_n
        print(n-sim_undone)

        z=np.random.normal(0, 1, (n_genes, cur_n))
        #print(z.shape)
        mzscores_tile = np.transpose(np.tile(mean_zscores, (cur_n, 1)))
        #print(mzscores_tile.shape)
        sim_zscores = mzscores_tile + np.dot(e_matrix, z)
        #print(sim_zscores.shape)
        sim_pvalues = np.transpose(norm.cdf(sim_zscores))
        #print(sim_pvalues.shape)
        for sim in sim_pvalues:
            #print(len(sim))
            cpma = calculate_cpma(sim, n_genes)
            sim_cpma.append(cpma)

    print('simulated cpma calculated')

    sim_cpma = np.array(sim_cpma)
   # np.savetxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_simzscores100.csv', sim_zscores, delimiter='\t')
    np.savetxt(output, sim_cpma, delimiter='\t', fmt='%f')


def run_decomposition_sim(in_zscores, out_sim, nsim):
    t0 = time.time()
    cov, mzscores = getCov(in_zscores)
    t1 = time.time()
    print(f'Cov and mean zscores calculation: {t1-t0}')
    e_values, Q = calculate_values(cov)
    t2 = time.time()
    print(f'Eigendecomposition: {t2-t1}')
    cov = 0
    simulate_cpma(mzscores, e_values, Q, out_sim, nsim)
    t3 = time.time()
    print(f'CPMA simulations: {t3-t2}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input file with zscores from matrix eQTL")
    parser.add_argument("-o", "--output", required=True, help="Ouptput file with simulated cpma values")
    parser.add_argument("-n", "--simulations", required=True, type=int, help="Number of simulated cpma values")
    params = parser.parse_args()
    run_decomposition_sim(params.input, params.output, params.simulations)


if __name__ == "__main__":
    main()
