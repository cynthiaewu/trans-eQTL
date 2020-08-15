import numpy as np
from scipy import linalg as LA
import argparse

def calculate_values(covfile, e_out, q_out):
    #zscores = genfromtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_zscores.csv', delimiter='\t', skip_header=1)
    #cov = np.loadtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_cov.csv', delimiter='\t', skiprows=1)
    cov = np.loadtxt(covfile, delimiter='\t', skiprows=1)
    print('files read')

    e_values, Q = LA.eig(cov)
    e_values = e_values.real
    e_values[e_values < 0] = 0
    Q = Q.real
    #np.savetxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_evalues.csv', e_values, delimiter='\t')
    #np.savetxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_Q.csv', Q, delimiter='\t')
    np.savetxt(e_out, e_values, delimiter='\t')
    np.savetxt(q_out, Q, delimiter='\t')
    #diag_e_values = np.diag(e_values)
    #E = np.sqrt(diag_e_values)
    print('eigendecomposition values calculated')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--cov", required=True, help="Input covariance matrix file")
    parser.add_argument("-e", "--eigenvalues", required=True, help="Ouptput file with eigenvalues")
    parser.add_argument("-q", "--eigenvectors", required=True, help="Ouptput file with eigenvectors")
    params = parser.parse_args()
    calculate_values(params.cov, params.eigenvalues, params.eigenvectors)


if __name__ == "__main__":
    main()
