import numpy as np
import pickle
from scipy.sparse import csr_matrix, save_npz
import sys
import os

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


proteins = set()
for root, _, fileNames in os.walk(path + 'features/'):
	pass

if excludeUnannotated:
	fileNames = [f for f in fileNames if '+unannotated' not in f]
else:
	fileNames = [f for f in fileNames if '+unannotated' in f]

for filename in fileNames:

	with open(root + filename) as f:
		for line in f:
			fields = line.split()
			proteins.add(fields[0])
			proteins.add(fields[1])



with open('../data/' + species + '/annotations/' + prefix + 'geneNames.pkl', 'rb') as f:
	geneNamesY = pickle.load(f)

if excludeUnannotated:
	protein2row = dict()
	for i, g in enumerate(geneNamesY):
		protein2row[g] = i
else:

	with open('../data/' + species + '/interactions/unannotatedProteinOrderBiogrid.pkl', 'rb') as f:
		protein2row = pickle.load(f)

	i = np.max([protein2row[k] for k in protein2row])
	for p in proteins:
		if p not in protein2row:
			i += 1
			protein2row[p] = i


for filename in fileNames:

	A = np.zeros((len(protein2row), len(protein2row)), float)

	with open(root + filename) as f:
		for line in f:
			[p1, p2, value] = line.split()


			i = protein2row[p1]
			j = protein2row[p2]
			A[i, j] = float(value)
			A[j, i] = float(value)


	if np.max(A) > 0.:

		save_npz(path + 'final/string/A_' + filename + '.npz', csr_matrix(A))

	else:
		print(filename)
