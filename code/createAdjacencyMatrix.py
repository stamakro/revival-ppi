import numpy as np
import pickle
from scipy.sparse import csr_matrix, save_npz
import sys

species = sys.argv[1]
try:
	prefix = sys.argv[2]
except IndexError:
	prefix = 'P_'

path = '../data/' + species + '/interactions/'

proteins = set()
with open(path + 'ppi-clean') as f:
	for line in f:
		for protein in line.split():
			proteins.add(protein)

with open('../data/' + species + '/annotations/' + prefix + 'geneNames.pkl', 'rb') as f:
	geneNamesY = pickle.load(f)

protein2row = dict()
for i, g in enumerate(geneNamesY):
	protein2row[g] = i


A = np.zeros((len(protein2row), len(protein2row)), int)


with open(path + 'ppi-clean') as f:
	for line in f:
		[p1, p2] = line.split()

		i = protein2row[p1]
		j = protein2row[p2]
		A[i, j] = 1
		A[j, i] = 1



save_npz(path + 'final/biogrid/A.npz', csr_matrix(A))
