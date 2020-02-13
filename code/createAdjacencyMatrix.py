import numpy as np
#import matplotlib.pyplot as plt
import pickle
from scipy.sparse import csr_matrix, save_npz
import sys

species = sys.argv[1]
path = '../data/' + species + '/interactions/'


proteins = set()
with open(path + 'ppi-clean') as f:
	for line in f:
		for protein in line.split():
			proteins.add(protein)


A = np.zeros((len(proteins), len(proteins)), int)

protein2row = dict()
row2protein = dict()

currentMax = 0

with open(path + 'ppi-clean') as f:
	for line in f:
		[p1, p2] = line.split()


		if p1 not in protein2row:
			protein2row[p1] = currentMax
			row2protein[currentMax] = p1
			currentMax += 1


		if p2 not in protein2row:
			protein2row[p2] = currentMax
			row2protein[currentMax] = p2
			currentMax += 1

		i = protein2row[p1]
		j = protein2row[p2]
		A[i, j] = 1
		A[j, i] = 1



save_npz(path + 'final/biogrid/A.npz', csr_matrix(A))

with open(path + 'final/biogrid/row2protein.pkl','wb') as f:
	pickle.dump(row2protein, f)



#sanity check - scale-free distribution
'''
degrees = np.sum(A, axis=0)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(degrees, bins=30,edgecolor='k')

plt.show()
'''
