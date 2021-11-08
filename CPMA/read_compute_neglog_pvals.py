import gzip
from collections import defaultdict
import numpy as np
from scipy import stats
import json
import argparse


def calculate_neglogpvals(input, output):
    snps = defaultdict(int)
    #with gzip.open('/storage/cynthiawu/trans_eQTL/GTex.v8/Run4_080221/Nerve-Tibial/chr21/matrixeqtl_results_ciseqtls_chr21.gz','r') as f:   
    print(f'Starting {input}')
    with gzip.open(input,'r') as f:        
        header = f.readline().decode("utf-8")
        for line in f:
            line = line.decode("utf-8")
            values = line.strip().split('\t')
            pval = 2*stats.norm.cdf(-np.abs(float(values[3])))
            snps[values[0]] += -np.log(pval)
    print(f'Finished dictionary of {input}')
    #with open('/storage/cynthiawu/trans_eQTL/GTex.v8/Run4_080221/Nerve-Tibial/chr21/neglogpvals.json', 'w') as outfile:
    with open(output, 'w') as outfile:
        json.dump(snps, outfile)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input matrix eqtl gz file")
    parser.add_argument("-o", "--output", required=True, help="Ouptput file with -log(pvals) for all snps")
    params = parser.parse_args()
    calculate_neglogpvals(params.input, params.output)


if __name__ == "__main__":
    main()
