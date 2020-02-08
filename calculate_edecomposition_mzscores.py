import numpy as np
from numpy import linalg as LA
from numpy import genfromtxt

def calculate_values(zfile, covfile, m_out, e_out, q_out):
    #zscores = genfromtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_zscores.csv', delimiter='\t', skip_header=1)
    #cov = np.loadtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_cov.csv', delimiter='\t', skiprows=1)
    zscores = genfromtxt(zfile, delimiter='\t', skip_header=1)
    cov = np.loadtxt(covfile, delimiter='\t', skiprows=1)
    print('files read')

    zscores_trans = np.transpose(zscores)
    num_genes = len(zscores_trans)
    mean_zscores = []
    for i in range(1,num_genes):
        mean_zscores.append(np.mean(zscores_trans[i][1:]))
    mean_zscores = np.array(mean_zscores)
    #np.savetxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_meanzscores.csv', mean_zscores, delimiter='\t')
    np.savetxt(m_out, mean_zscores, delimiter='\t')
    print('mean zscores calculated')

    e_values, Q = LA.eig(cov)
    #np.savetxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_evalues.csv', e_values, delimiter='\t')
    #np.savetxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_Q.csv', Q, delimiter='\t')
    np.savetxt(e_out, e_values, delimiter='\t')
    np.savetxt(q_out, Q, delimiter='\t')
    #diag_e_values = np.diag(e_values)
    #E = np.sqrt(diag_e_values)
    print('eigendecomposition values calculated')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-z", "--zscores", required=True, help="Input zscores file")
    parser.add_argument("-c", "--cov", required=True, help="Input covariance matrix file")
    parser.add_argument("-m", "--mzscores", required=True, help="Ouptput file with mean zscores")
    parser.add_argument("-e", "--eigenvalues", required=True, help="Ouptput file with eigenvalues")
    parser.add_argument("-q", "--eigenvectors", required=True, help="Ouptput file with eigenvectors")
    params = parser.parse_args()
    calculate_values(params.zscores, params.cov, params.mzscores, params.eigenvalues, params.eigenvectors)


if __name__ == "__main__":
    main()
