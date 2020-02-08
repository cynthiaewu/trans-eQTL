import pandas as pd
import numpy as np
from numpy import genfromtxt

def calculate_cpma(input, output):
    #pvalues = genfromtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_pvalues_nofdr.csv', delimiter='\t', skip_header=1)
    pvalues = genfromtxt(input, delimiter='\t', skip_header=1)
    #np.negative(np.log(pvalues[1:]))
    num_snps = len(pvalues)
    num_genes = len(pvalues[1])-1
    cpma = []
    # 3623 snps
    for i in range(num_snps):
        likelihood = np.mean(np.negative(np.log(pvalues[i][1:])))
        value = -2 * ((((likelihood - 1) * num_genes)/likelihood) - num_genes*np.log(likelihood))
        cpma.append(value)

    snps = pd.read_csv(input, usecols=[0], sep='\t', names = ['snp'], skiprows = 1)
    snps.insert(column = 'cpma', loc = 1, value = cpma)
    #snps.to_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_cpma_nofdr.csv', index=False, header=True, sep='\t')
    snps.to_csv(output, index=False, header=True, sep='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input pvalues no fdr file")
    parser.add_argument("-o", "--output", required=True, help="Ouptput file with cpma values")
    params = parser.parse_args()
    calculate_cpma(params.input, params.output)


if __name__ == "__main__":
    main()
