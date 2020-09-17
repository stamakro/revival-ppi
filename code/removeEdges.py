import numpy as np
import pickle
import sys
from scipy.spatial.distance import cdist
from scipy.sparse import csr_matrix
from routinesnetworks import *
from routinesgo import semanticDistance
import os
from copy import deepcopy

def removeEdges(A, fractionToRemove, mode='random', seed=104):
    edges = np.where(A)
    edgePairs = np.array([(i,j) for i, j in zip(edges[0], edges[1]) if i < j])

    np.random.seed(seed)
    toremove = int(np.round(len(edgePairs) * fractionToRemove))
    if mode == 'random':   

        pairsToRemove = np.random.permutation(edgePairs)[:toremove]



    elif mode == 'degree':
        degrees = np.sum(A,0)
        edgeDeg = np.zeros(edgePairs.shape[0])
        
        for i, (n1, n2) in enumerate(edgePairs):
            edgeDeg[i] = np.minimum(degrees[n1], degrees[n2])

        keepProba = edgeDeg / np.sum(edgeDeg)

        indToKeep = np.random.choice(np.arange(edgePairs.shape[0]), size=edgePairs.shape[0] - toremove, replace=False, p=keepProba)             
        indToRemove = np.setdiff1d(np.arange(edgePairs.shape[0]), indToKeep)
        pairsToRemove = edgePairs[indToRemove]

    else:
        print('Use mode=\'random\' or \'degree\'')
        sys.exit(1)


    Anew = deepcopy(A)
    for p in pairsToRemove:
        Anew[p[0], p[1]] = 0
        Anew[p[1], p[0]] = 0


    return Anew


experimentPath = '../experiments/yeast_P_median/'
stuff = experimentPath.split('/')[2]
stuff2 = stuff.split('_')
species = stuff2[0]

annoPrefix = stuff2[1] + '_'



print('Loading go...')
[Y, geneNames, termNames, gene2row, term2col] = goLoader(species)

with open(experimentPath + 'icDict.pkl', 'rb') as f:
    icDict = pickle.load(f)

ic = np.array([icDict[t] for t in termNames])

print('Loading experimental PPI network...')

Aexp = getPPInetwork(species, 'biogrid')
assert np.max(np.abs((Aexp - Aexp.T))) < 1e-10

try:
    Aexp2 = getPPInetwork(species, 'experiments')
    assert np.max(np.abs((Aexp2 - Aexp2.T))) < 1e-10

    if np.max(Aexp2) > 0.:
    	threshold = np.median(Aexp2[Aexp2 > 0])
    	Aexp2 = (Aexp2 >= threshold).astype(int)
    	Aexp = np.maximum(Aexp, Aexp2)
except FileNotFoundError:
    pass

Aexp = Aexp.astype(float)

seeds = [121, 898921630, 283261, 2313,51437809]
fractions = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99])

try:
    os.mkdir('../data/yeast/networks/downsampled/')
except FileExistsError:
    pass

for i, s in enumerate(seeds):

    try:
        os.mkdir('../data/yeast/networks/downsampled/seed' + str(i))
    except FileExistsError:
        pass


    for frac in fractions:
        print('%d, %f' % (i, frac))
        Ar = removeEdges(Aexp, frac, 'random', s)
        Ad = removeEdges(Aexp, frac, 'degree', s)

        with open('../data/yeast/networks/downsampled/seed' + str(i) + '/random_' + str(frac) + '.pkl', 'wb') as f:
            pickle.dump(csr_matrix(Ar), f)

        ee = np.where(Ar)
        edges = [(i,j) for (i,j) in zip(ee[0], ee[1]) if i < j]

        with open('../data/yeast/networks/downsampled/seed' + str(i) + '/random_' + str(frac) + '.network', 'w') as f:
            for e in edges:
                f.write(str(e[0]) + '\t' + str(e[1]) + '\n')

        with open('../data/yeast/networks/downsampled/seed' + str(i) + '/degree_' + str(frac) + '.pkl', 'wb') as f:
            pickle.dump(csr_matrix(Ad), f)

        ee = np.where(Ad)
        edges = [(i,j) for (i,j) in zip(ee[0], ee[1]) if i < j]

        with open('../data/yeast/networks/downsampled/seed' + str(i) + '/degree_' + str(frac) + '.network', 'w') as f:
            for e in edges:
                f.write(str(e[0]) + '\t' + str(e[1]) + '\n')

