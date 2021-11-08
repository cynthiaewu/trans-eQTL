import gzip
from collections import defaultdict
import numpy as np
from scipy import stats
import json
import argparse


def calculate_cpma(input, output):
    #with open('/storage/cynthiawu/trans_eQTL/GTex.v8/Run4_080221/Nerve-Tibial/chr21/neglogpvals.json') as json_file:
    with open(input) as json_file:
        data = json.load(json_file)

    num_genes = 17156
    cpma = defaultdict(int)
    for key, value in data.items():
        if value == 0:
            cpma[key] = 0
            #print(key, ((value)/num_genes))
        else:
            likelihood = 1/((value)/num_genes)
            cpma_value = -2 * ((((likelihood - 1) * num_genes)/likelihood) - num_genes*np.log(likelihood))
            cpma[key] = cpma_value
    #with open('/storage/cynthiawu/trans_eQTL/GTex.v8/Run4_080221/Nerve-Tibial/chr21/cpma_chr21', 'w') as f:
    with open(output, 'w') as f:
        f.write(f'SNP\tCPMA\n')
        for key, value in dict(sorted(cpma.items(), reverse=True, key=lambda item: item[1])).items():
            f.write(f'{key}\t{value}\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input -log(pvals) for all snps")
    parser.add_argument("-o", "--output", required=True, help="Output cpma file")
    params = parser.parse_args()
    calculate_cpma(params.input, params.output)


if __name__ == "__main__":
    main()
