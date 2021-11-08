import numpy as np
import pandas as pd
from scipy import linalg as LA
import os
import argparse

def getCov(input, m_out):
    print(f'Reading datafile {input}')
    data = pd.read_csv(input, sep='\t', index_col=0)
    genes = list(data.columns)
    scores = np.transpose(np.array(data))

    cov = np.cov(scores)
    cov_out = '/storage/cynthiawu/trans_eQTL/GTex.v8/Run3_070221/Adipose-Subcutaneous_Euro/Adipose-Subcutaneous_cov.txt'
    if not os.path.isfile(cov_out):
        np.savetxt(cov_out, cov, delimiter='\t', fmt='%f')
    #cov_matrix = pd.DataFrame(cov, columns=genes)
    print('cov matrix calculated')

    mean_zscores = []
    for i in scores:
        mean_zscores.append(np.mean(i))
    mean_zscores = np.array(mean_zscores)
    np.savetxt(m_out, mean_zscores, delimiter='\t', fmt='%f')
    print('mean zscores calculated')

    return cov


def calculate_values(input, m_out, e_out, q_out):
    #zscores = genfromtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_zscores.csv', delimiter='\t', skip_header=1)
    #cov = np.loadtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_cov.csv', delimiter='\t', skiprows=1)
    #cov = np.loadtxt(covfile, delimiter='\t', skiprows=1)
    cov = getCov(input, m_out)
    #print(f'files read {covfile}')

    e_values, Q = LA.eigh(cov)
    e_values = e_values.real
    e_values[e_values < 0] = 0
    Q = Q.real
    #np.savetxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_evalues.csv', e_values, delimiter='\t')
    #np.savetxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_Q.csv', Q, delimiter='\t')
    np.savetxt(e_out, e_values, delimiter='\t')
    np.savetxt(q_out, Q, delimiter='\t')
    #diag_e_values = np.diag(e_values)
    #E = np.sqrt(diag_e_values)
    print(f'eigendecomposition values calculated for {input}')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input file with zscores from matrix eQTL")
    parser.add_argument("-m", "--m_out", required=True, help="Ouptput file with mean zscores for genes")
    parser.add_argument("-e", "--eigenvalues", required=True, help="Ouptput file with eigenvalues")
    parser.add_argument("-q", "--eigenvectors", required=True, help="Ouptput file with eigenvectors")
    params = parser.parse_args()
    calculate_values(params.input, params.m_out,  params.eigenvalues, params.eigenvectors)


if __name__ == "__main__":
    main()
