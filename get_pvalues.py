import pandas as pd
import numpy as np

data = pd.read_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls.csv', sep='\t')
data_sorted = data.sort_values(["SNP", "gene"], ascending=[True, True])
genes = list(data_sorted[0:28315]['gene'])

cur_snp = ''
cur_row = []
scores = []
snps = []
count = 0
for index, row in data_sorted.iterrows():
    new_snp = row[0]
    pvalue = row[4]
    if cur_snp == new_snp:
        cur_row.append(pvalue)
    else:
        scores.append(cur_row)
        cur_row = [pvalue]
        cur_snp = new_snp
        snps.append(new_snp)
        #print(new_snp)
        count += 1
        if count % 100 == 0:
            print(count)
scores.append(cur_row)
scores.remove(scores[0])

pvalue_matrix = pd.DataFrame(scores, index=snps, columns=genes)
pvalue_matrix.to_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_pvalues_nofdr.csv', index=True, header=True, sep='\t')

