import pandas as pd
import numpy as np

data = pd.read_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls.csv', sep='\t')
data_sorted = data.sort_values(["SNP", "gene"], ascending=[True, True])
genes = list(data_sorted[0:28315]['gene'])

cur_snp = ''
cur_row = []
scores = []
snps = []
#start
count = 0
for index, row in data_sorted.iterrows():
    new_snp = row[0]
    zscore = row[3]
    if cur_snp == new_snp:
        cur_row.append(zscore)
    else:
        scores.append(cur_row)
        cur_row = [zscore]
        cur_snp = new_snp
        snps.append(new_snp)
        #print(new_snp)
        count += 1
        if count % 100:
            print(count)
scores.append(cur_row)
scores.remove(scores[0])

zscore_matrix = pd.DataFrame(scores, index=snps, columns=genes)
zscore_matrix.to_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_zscores.csv', index=True, header=True, sep='\t')

cov = np.cov(np.transpose(scores))
cov_matrix = pd.DataFrame(cov, columns=genes)
cov_matrix.to_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_cov.csv', index=True, header=True, sep='\t')
