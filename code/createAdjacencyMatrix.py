import numpy as np
import pickle
from scipy.sparse import csr_matrix, save_npz
import sys

species = sys.argv[1]
try:
	prefix = sys.argv[2]
except IndexError:
	prefix = 'P_'

try:
	excludeUnannotated = bool(int(sys.argv[3]))
except IndexError:
	excludeUnannotated = True

path = '../data/' + species + '/interactions/'

if excludeUnannotated:
	inFileName = path + 'ppi-clean'
else:
	inFileName = path + 'ppi-clean+unannotated'


proteins = set()
with open(inFileName) as f:
	for line in f:
		for protein in line.split():
			proteins.add(protein)

with open('../data/' + species + '/annotations/' + prefix + 'geneNames.pkl', 'rb') as f:
	geneNamesY = pickle.load(f)

protein2row = dict()
for i, g in enumerate(geneNamesY):
	protein2row[g] = i

if not excludeUnannotated:
	for p in proteins:
		if p not in protein2row:
			i += 1
			protein2row[p] = i

	with open('../data/' + species + '/interactions/unannotatedProteinOrderBiogrid.pkl', 'wb') as f:
		pickle.dump(protein2row, f)

A = np.zeros((len(protein2row), len(protein2row)), int)


with open(inFileName) as f:
	for line in f:
		[p1, p2] = line.split()

		i = protein2row[p1]
		j = protein2row[p2]
		A[i, j] = 1
		A[j, i] = 1


outFile = path + ''


if excludeUnannotated:
	save_npz(path + 'final/biogrid/A.npz', csr_matrix(A))

else:
	save_npz(path + 'final/biogrid/A+unannotated.npz', csr_matrix(A))
