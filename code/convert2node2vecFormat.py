from scipy.sparse import load_npz
from routinesnetworks import getPPInetwork
import pickle
import sys
import numpy as np
from copy import deepcopy

species = sys.argv[1]

path = '../data/' + species + '/interactions/final/biogrid/'

A = load_npz(path + 'A+unannotated.npz').toarray()

Aexp2 = getPPInetwork(species, 'experiments+unannotated')

if np.max(Aexp2) > 0.:

    threshold = np.median(Aexp2[Aexp2 > 0])

    Aexp2 = (Aexp2 >= threshold).astype(int)

    Afinal = deepcopy(Aexp2)

    Afinal[:A.shape[0], :A.shape[0]] = np.maximum(A, Aexp2[:A.shape[0], :A.shape[0]])

else:
    Afinal = A


ii, jj = np.where(Afinal)

with open(path + 'netIncludingUnannotated.network', 'w') as f:
	for ind1, ind2 in zip(ii, jj):
		if ind1 > ind2:
			f.write(str(ind1) + '\t' + str(ind2) + '\n')
