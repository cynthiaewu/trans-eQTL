import numpy as np
from numpy import linalg as LA
from numpy import genfromtxt

mean_zscores = np.loadtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_meanzscores.csv', dtype=complex, delimiter='\t')
print('mean zscores file read')

e_values = np.loadtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_evalues.csv', dtype=complex, delimiter='\t')
print('e_values file read')
Q = np.loadtxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_Q.csv', dtype=complex, delimiter='\t')
print('Q file read')
diag_e_values = np.diag(e_values)
E = np.sqrt(diag_e_values)

print('starting simulations')
sim_zscores = []
e_matrix = np.dot(Q, E)
for i in range(100):
    if (i%10==0):
        print(i)
    z = np.random.normal(0, 1, 28315)
    cur_sim_zscores = mean_zscores + np.dot(e_matrix,z)
    sim_zscores.append(cur_sim_zscores.real)
print('simuated zscores calculated')

sim_zscores = np.array(sim_zscores)
np.savetxt('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr1_gene_snp_eqtls_simzscores100.csv', sim_zscores, delimiter='\t')
