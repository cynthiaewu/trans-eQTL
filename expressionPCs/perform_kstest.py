import pandas as pd
import numpy as np
from scipy import stats
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


def iter_simulations(input, method, topx, iterations, num_test):
    sim_prefix = 'Simulation'
    for i in range(iterations):
        if method==0:
            in_file = f'{input}/{sim_prefix}_{i}/CPMA/gene-snp-eqtl_cpma_fixed'
            out_file = f'{input}/{sim_prefix}_{i}/CPMA/gene-snp-eqtl_cpma_pvalues_fixed'
        if method==1:
            in_file = f'{input}/{sim_prefix}_{i}/CPMAx/gene-snp-eqtl_cpma_topx_{topx}'
            out_file = f'{input}/{sim_prefix}_{i}/CPMAx/gene-snp-eqtl_cpmax_pvalues_{topx}'
        if method==2:
            in_file = f'{input}/{sim_prefix}_{i}/expressionPCs/gene-snp-eqtl_PCs_cpma_topx_{topx}'
            out_file = f'{input}/{sim_prefix}_{i}/expressionPCs/gene-snp-eqtl_PCs_cpmax_pvalues_{topx}'
        compute_pvalue(in_file, num_test, out_file)
    print(f'finished for {input}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input file with cpma values per snp")
    parser.add_argument("-m", "--method", type=int, required=True, help="0 for cpma, 1 for cpma_topx, or 2 for cpma_topx_PCs")
    parser.add_argument("-x", "--topx", default=0.1, type=float, help="Top x percent of genes to be used for cpma calculation")
    parser.add_argument("-n", "--iterations", type=int, required=True, help="# iterations to simulate genotype and expression files")
    parser.add_argument("-t", "--num_test", type=int, required=True, help="Number of tests to correct pvalues with")
    #parser.add_argument("-o", "--output", required=True, help="Output file with cpma values, pvalues, adjusted pvalues per snp")
    params = parser.parse_args()

    iter_simulations(input=params.input,
          method=params.method,
          topx=params.topx,
          iterations=params.iterations,
          num_test=params.num_test)
          #output=params.output)


if __name__ == "__main__":
    main()
