import pandas as pd
from scipy.stats import chi2
import argparse


#compute pvalue from chi distribution with df=1
def compute_pvalue(input, num_test, output):
    df = pd.read_csv(input, sep='\t', index_col=0)
    pvalue = []
    adj_pvalue = []
    for index, row in df.iterrows():
        cpma = row['cpma']
        value = 1 - chi2.cdf(cpma, 1)
        pvalue.append(value)
        adj_pvalue.append(min(1, value * num_test))
    df['pvalue'] = pvalue
    df['adj_pvalue'] = adj_pvalue
    df.to_csv(output, sep='\t')


def iter_simulations(input, iterations, num_test):
    sim_prefix = 'Simulation'
    for i in range(iterations):
        in_file = f'{input}/{sim_prefix}_{i}/CPMA/gene-snp-eqtl_cpma'
        out_file = f'{input}/{sim_prefix}_{i}/CPMA/gene-snp-eqtl_cpma_pvalues'
        compute_pvalue(in_file, num_test, out_file)
    print(f'finished for {input}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input file with cpma values per snp")
    parser.add_argument("-n", "--iterations", type=int, required=True, help="# iterations to simulate genotype and expression files")
    parser.add_argument("-t", "--num_test", type=int, required=True, help="Number of tests to correct pvalues with")
    #parser.add_argument("-o", "--output", required=True, help="Output file with cpma values, pvalues, adjusted pvalues per snp")
    params = parser.parse_args()

    iter_simulations(input=params.input,
          iterations=params.iterations,
          num_test=params.num_test)
          #output=params.output)


if __name__ == "__main__":
    main()
