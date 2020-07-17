import argparse
import numpy as np
from scipy.stats import norm


def calculate_cpma(sim_pvalues, num_genes):
    likelihood = np.mean(np.negative(np.log(sim_pvalues)))
    value = -2 * ((((likelihood - 1) * num_genes)/likelihood) - num_genes*np.log(likelihood))
    return value

def get_cpma(simfile, output):
    sim_zscores = np.loadtxt(simfile, delimiter='\t') 
    n_genes = len(sim_zscores[0])
    print(n_genes)
    sim_cpma = []
    sim_pvalues = norm.cdf(sim_zscores)
    #print(sim_pvalues.shape)
    for sim in sim_pvalues:
        cpma = calculate_cpma(sim, n_genes)
        sim_cpma.append(cpma)

    sim_cpma = np.array(sim_cpma)
   # np.savetxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_simzscores100.csv', sim_zscores, delimiter='\t')
    np.savetxt(output, sim_cpma, delimiter='\t', fmt='%f')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-z", "--simzscores", required=True, help="Input sim zscores file")
    parser.add_argument("-o", "--output", required=True, help="Ouptput file with simulated cpma values")
    params = parser.parse_args()
    get_cpma(params.simzscores, params.output)


if __name__ == "__main__":
    main()
