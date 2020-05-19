from routinesnetworks import getPPInetwork
import sys
import numpy as np


species = sys.argv[1]

Aexp = getPPInetwork(species, 'biogrid')
assert np.max(np.abs((Aexp - Aexp.T)) < 1e-10)

try:
    Aexp2 = getPPInetwork(species, 'experiments')
    assert np.max(np.abs((Aexp2 - Aexp2.T))) < 1e-10

    if np.max(Aexp2) > 0.:
    	threshold = np.median(Aexp2[Aexp2 > 0])
    	Aexp2 = (Aexp2 >= threshold).astype(int)
    	Aexp_total = np.maximum(Aexp, Aexp2)
except FileNotFoundError:
    Aexp2 = None
    Aexp_total = Aexp

Aexp_total = Aexp_total.astype(float)

N = Aexp_total.shape[0]

print(N)

Nedges = np.sum(Aexp_total) / 2

print(Nedges)

Npairs = N * (N-1) / 2

print(Npairs)

print(100 * Nedges / Npairs)

empty = np.mean(np.sum(Aexp_total,0) == 0) * 100

print(empty)



if species == 'tomato':
	datasources = np.array(['neighborhood_transferred', 'coexpression_transferred', 'experiments_transferred', 'textmining', 'cooccurence','fusion', 'textmining_transferred', 'homology'])

elif species == 'ecoli':
	datasources = np.array(['neighborhood', 'neighborhood_transferred', 'coexpression_transferred', 'experiments_transferred', 'textmining', 'cooccurence','fusion', 'textmining_transferred', 'homology'])


else:
	datasources = np.array(['coexpression', 'neighborhood_transferred', 'coexpression_transferred', 'experiments_transferred', 'textmining', 'cooccurence','fusion', 'textmining_transferred', 'homology'])


AA = np.zeros((2+datasources.shape[0], Aexp.shape[0], Aexp.shape[1]))


for i, ds in enumerate(datasources,2):
	print ('Reading predicted %d' % i)
	At = getPPInetwork(species, ds)
	AA[i] = At 

AA[0] = Aexp

if Aexp2 is not None:
	AA[1] = Aexp2


print('Correcting for homology...')

#homology correction
homInd = np.where(datasources == 'homology')[0] + 2

tmInd = np.where(datasources == 'textmining')[0] + 2

AA[tmInd] = AA[tmInd] * (1 - AA[homInd])

coInd = np.where(datasources == 'cooccurence')[0] + 2
AA[coInd] = AA[coInd] * (1 - AA[homInd])

for i in range(2, AA.shape[0]):
	At = AA[i]

	thres = np.median(At[At > 0])

	AA[i] = (At >= thres).astype(int)



datasourcesAll = ['biogrid', 'experiments'] + list(datasources) 


for i in range(5):
	print('')


tri, trj = np.triu_indices_from(AA[0], k=1)
allInd = np.arange(len(datasourcesAll))

for i,ds in enumerate(datasourcesAll):
	currentInd = np.delete(allInd, i)
	print(ds)
	
	currentNet = AA[i]

	restCombined = np.max(AA[currentInd], axis=0)

	ee = np.where(currentNet[tri, trj])[0]
	

	print('Total edges: %d'  % ee.shape[0])

	restFlat = restCombined[tri, trj]



	print('Unique edges: %d' % np.sum(restFlat[ee] == 0))
	print('')





