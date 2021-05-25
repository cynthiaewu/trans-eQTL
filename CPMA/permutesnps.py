import pandas as pd
import numpy as np
import random
import argparse


def permute_snps(input, snp, n, output):
#     snp_geno = np.array(genotype.loc['SNP0'])
    genotype = pd.read_csv(input, sep='\t', index_col=0)
    snp_geno = np.array(genotype.loc[snp])
    permute_snp = []
    for i in range(n):
        permute_snp.append(np.random.permutation(snp_geno))
    row = ['SNPermute' + str(i) for i in range(n)]
    permute_snp_df = pd.DataFrame(permute_snp, index=row, columns=genotype.columns)
    permute_snp_df.to_csv(output, sep='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input genotype file")
    parser.add_argument("-s", "--snp", required=True, help="Name of snp to permute")
    parser.add_argument("-n", "--n_permutations", required=True, type=int, help="Number of times to permute snp")
    parser.add_argument("-o", "--output", required=True, help="Output file permuted genotypes of specified snp")
    params = parser.parse_args()
    permute_snps(params.input, params.snp, params.n_permutations, params.output)


if __name__ == "__main__":
    main()
