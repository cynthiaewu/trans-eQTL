import pandas as pd
import sklearn.decomposition
import gzip
import argparse


def get_sigresults(folder, method, factor):
    #sig_results = open('/storage/cynthiawu/trans_eQTL/GTex.v8/Run2_031621/Nerve-Tibial_Euro/all_chr/matrixeQTL_PCs_sig_results', 'a')
    
    sig_results = open(f'{folder}/matrixeQTL_PCs_sig_results_cpma', 'a')
    if method == 1:
        sig_results = open(f'{folder}/matrixeQTL_PCs_sig_results_pairwise_{factor}factors', 'at')
    threshold = 0.05
    for i in range(1, 23):
        #path = f'{folder}/chr{i}/expressionPCs/matrixeQTL_PCs_cpma_converted_pval_chr{i}'
        path = f'{folder}/chr{i}/CPMA/gene-snp_eqtls_empiricalpvalues_chr{i}_converted'
        if method == 1:
        #path = f'{folder}/chr{i}/expressionPCs/matrixeQTL_PCs_chr{i}'
            path = f'{folder}/chr{i}/matrixeqtl/matrixeQTL_results_chr{i}_{factor}factors.gz'
        #with open(path) as f:
        with gzip.open(path, 'rb') as f:
            header = f.readline().decode("utf-8") 
            #header = f.readline() 
            
            # write header only once, from chr1
            if i == 1:
                sig_results.write(header)
            for line in f:
                line = line.decode("utf-8") 
                values = (line.strip().split('\t'))
                if method == 1:
                    if float(values[4]) > threshold:
                       break
                    sig_results.write(line)
                if method == 2:
                    if float(values[2]) <= threshold: 
                        sig_results.write(line)
            
 
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", required=True, help="Input folder with all chr folders")
    parser.add_argument("-m", "--method", required=True, type=int, help="Method=1 for matrixeqtl pairwise results, Method=2 for cpma results")
    parser.add_argument("-k", "--factor",  type=int, help="# peer factors")
    params = parser.parse_args()
    get_sigresults(params.folder, params.method, params.factor)


if __name__ == "__main__":
    main()
