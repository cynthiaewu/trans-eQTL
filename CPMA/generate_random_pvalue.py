import random
import numpy as np

ran_matrix = np.random.randint(1,1001, size=(1000, 1000))
ran_matrix = ran_matrix*0.001
#print(ran_matrix)

np.savetxt('/storage/cynthiawu/trans_eQTL/test/ran_pvalue_matrix.csv', ran_matrix, fmt='%1.3f', delimiter='\t')
