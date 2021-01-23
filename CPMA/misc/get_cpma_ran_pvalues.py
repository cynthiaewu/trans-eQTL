#import pandas as pd
import numpy as np
from numpy import genfromtxt

pvalues = genfromtxt('/storage/cynthiawu/trans_eQTL/test/ran_pvalue_matrix.csv', delimiter='\t')
#np.negative(np.log(pvalues[1:]))
num_genes = 1000
cpma = []
for i in range(1000):
    likelihood = np.mean(np.negative(np.log(pvalues[i])))
    value = -2 * ((((likelihood - 1) * num_genes)/likelihood) - num_genes*np.log(likelihood))
    cpma.append(value)

output = '/storage/cynthiawu/trans_eQTL/test/ran_pvalue_cpma.csv'    
f = open(output, 'w')
for value in cpma:
    f.write('%f' % value + '\n')
f.close()

#snps = pd.read_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_pvalues.csv',usecols=[0], sep='\t', names = ['snp'], skiprows = 1)
#snps.insert(column = 'cpma', loc = 1, value = cpma)
#snps.to_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_cpma_nofdr.csv', index=False, header=True, sep='\t')
