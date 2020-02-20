import pandas as pd
import numpy as np
import argparse


def getCov(input, cov_out, m_out):
    data = pd.read_csv(input, sep='\t', index_col=0)
    genes = list(data.columns)
    scores = np.transpose(np.array(data))
    
    cov = np.cov(scores)
    cov_matrix = pd.DataFrame(cov, columns=genes)
    #cov_matrix.to_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_cov.csv', index=True, header=True, sep='\t')
    cov_matrix.to_csv(cov_out, index=False, header=True, sep='\t')
    print('cov matrix calculated')

    mean_zscores = []
    for i in scores:
        mean_zscores.append(np.mean(i))
    mean_zscores = np.array(mean_zscores)
    #np.savetxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_meanzscores.csv', mean_zscores, delimiter='\t')
    np.savetxt(m_out, mean_zscores, delimiter='\t', fmt='%f'))
    print('mean zscores calculated')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input file with zscores from matrix eQTL")
    parser.add_argument("-c", "--cov_out", required=True, help="Ouptput file with covariance matrix for genes")
    parser.add_argument("-m", "--m_out", required=True, help="Ouptput file with mean zscores for genes")
    params = parser.parse_args()
    getCov(params.input, params.cov_out, params.m_out)


if __name__ == "__main__":
    main()
