import pandas as pd
import numpy as np
import argparse

def getValues(input, num_genes, p_out, z_out): 
    #data = pd.read_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls.csv', sep='\t')
    data = pd.read_csv(input, sep='\t')
    data_sorted = data.sort_values(["SNP", "gene"], ascending=[True, True])
    genes = list(data_sorted[0:num_genes]['gene'])

    all_pvalues = []
    all_zscores = []

    snps = sorted(set(data_sorted['SNP']))

    total = len(data_sorted)
    iterations = int(total/num_genes)
    for i in range(iterations):
        current = i * num_genes
        all_pvalues.append(list(data_sorted[current:current+num_genes]['p-value']))
        #all_zscores.append(list(data_sorted[current:current+num_genes]['t-stat']))

    pvalue_matrix = pd.DataFrame(all_pvalues, index=snps, columns=genes)
    #pvalue_matrix.to_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_pvalues_nofdr.csv', index=True, header=True, sep='\t')
    pvalue_matrix.to_csv(p_out, index=True, header=True, sep='\t')
    #zscore_matrix = pd.DataFrame(all_zscores, index=snps, columns=genes)
    #zscore_matrix.to_csv(z_out, index=True, header=True, sep='\t')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input matrix eQTL file")
    parser.add_argument("-n", "--num_genes", type=int, required=True, help="Number of genes")
    parser.add_argument("-p", "--pvalues_out", required=True, help="Output file with matrix of pvalues for every snp and gene")
    parser.add_argument("-z", "--zscores_out", required=True, help="Output file with matrix of zscores for every snp and gene")
    params = parser.parse_args()
    getValues(params.input, params.num_genes, params.pvalues_out, params.zscores_out)


if __name__ == "__main__":
    main()
