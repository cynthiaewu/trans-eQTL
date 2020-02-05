import numpy as np
from numpy import linalg as LA
from numpy import genfromtxt

zscores = genfromtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_zscores.csv', delimiter='\t', skip_header=1)
cov = np.loadtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_cov.csv', delimiter='\t', skiprows=1)
print('files read')

zscores_trans = np.transpose(zscores)
mean_zscores = []
for i in range(1,28316):
    mean_zscores.append(np.mean(zscores_trans[i][1:]))
mean_zscores = np.array(mean_zscores)
np.savetxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_meanzscores.csv', mean_zscores, delimiter='\t')
print('mean zscores calculated')

e_values, Q = LA.eig(cov)
np.savetxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_evalues.csv', e_values, delimiter='\t')
np.savetxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_Q.csv', Q, delimiter='\t')
#diag_e_values = np.diag(e_values)
#E = np.sqrt(diag_e_values)
print('eigendecomposition values calculated')
