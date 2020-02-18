import pandas as pd
import numpy as np
import argparse

def getValues(input, p_out, z_out): 
    #data = pd.read_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls.csv', sep='\t')
    data = pd.read_csv(input, sep='\t')
    data_sorted = data.sort_values(["SNP", "gene"], ascending=[True, True])
    genes = list(data_sorted[0:28315]['gene'])

    cur_snp = ''
    cur_zrow = []
    cur_prow = []
    all_zscores = []
    all_pvalues = []
    snps = []
    count = 0
    for index, row in data_sorted.iterrows():
        new_snp = row[0]
        zscore = row[3]
        pvalue = row[4]
        if cur_snp == new_snp:
            cur_zrow.append(zscore)
            cur_prow.append(pvalue)
        else:
            all_pvalues.append(cur_prow)
            cur_prow = [pvalue]
            all_zscores.append(cur_zrow)
            cur_zrow = [zscore]
            cur_snp = new_snp
            snps.append(new_snp)
            #print(new_snp)
            count += 1
            if count % 1000 == 0:
                print(count)
    all_pvalues.append(cur_prow)
    all_pvalues.remove(all_pvalues[0])
    all_zscores.append(cur_zrow)
    all_zscores.remove(all_zscores[0])

    pvalue_matrix = pd.DataFrame(all_pvalues, index=snps, columns=genes)
    #pvalue_matrix.to_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_pvalues_nofdr.csv', index=True, header=True, sep='\t')
    pvalue_matrix.to_csv(p_out, index=True, header=True, sep='\t')
    zscore_matrix = pd.DataFrame(all_zscores, index=snps, columns=genes)
    zscore_matrix.to_csv(z_out, index=True, header=True, sep='\t')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input matrix eQTL file")
    parser.add_argument("-p", "--pvalues_out", required=True, help="Output file with matrix of pvalues for every snp and gene")
    parser.add_argument("-z", "--zscores_out", required=True, help="Output file with matrix of zscores for every snp and gene")
    params = parser.parse_args()
    getValues(params.input, params.pvalues_out, params.zscores_out)


if __name__ == "__main__":
    main()
