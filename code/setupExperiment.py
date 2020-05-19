from itertools import combinations
from scipy.special import binom
from sklearn.model_selection import KFold
import numpy as np
from routinesgo import *
from routinesnetworks import *
import os
import sys
from subprocess import call


species = sys.argv[1]

try:
    annoPrefix = sys.argv[2]
except IndexError:
    annoPrefix = 'P_'

#median or fixed
try:
    thresholdType = sys.argv[3]
except IndexError:
    thresholdType = 'median'

if thresholdType == 'fixed':
    try:
        threshold = float(sys.argv[4])
    except IndexError:
        threshold = 0.5

dag, mapping = read_ontology_from_file('../data/go/go-final.obo')

print('Loading go...')
[Y, geneNames, termNames, gene2row, term2col] = goLoader(species, annoPrefix)

print('Loading experimental PPI network...')

Aexp = getPPInetwork(species, 'biogrid')
try:
    Aexp2 = getPPInetwork(species, 'experiments')

    if np.max(Aexp2) > 0.:

        if thresholdType == 'median':
            threshold = np.median(Aexp2[Aexp2 > 0])

        Aexp2 = (Aexp2 >= threshold).astype(int)
        Aexp = np.maximum(Aexp, Aexp2)
except FileNotFoundError:
    pass

for root, _, names in os.walk('../data/' + species + '/interactions/final/string/'):

    datasources = []

    for n in names:
        fields = n.split('.')[0].split('_')[1:]
        if len(fields) == 1:
            if fields[0] != 'database' and fields[0] != 'experiments':
                datasources.append(fields[0])

        else:
            datasources.append(fields[0] + '_' + fields[1])

    datasources = np.array(datasources)



if species == 'tomato':
        datasources2 = np.array(['neighborhood_transferred', 'coexpression_transferred', 'experiments_transferred', 'textmining', 'cooccurence','fusion', 'textmining_transferred', 'homology'])

elif species == 'ecoli':
        datasources2 = np.array(['neighborhood', 'neighborhood_transferred', 'coexpression_transferred', 'experiments_transferred', 'textmining', 'cooccurence','fusion', 'textmining_transferred', 'homology'])


else:
        datasources2 = np.array(['coexpression', 'neighborhood_transferred', 'coexpression_transferred', 'experiments_transferred', 'textmining', 'cooccurence','fusion', 'textmining_transferred', 'homology'])

assert set(datasources) == set(datasources2)
datasources = datasources2




AA = np.zeros((datasources.shape[0], Aexp.shape[0], Aexp.shape[1]))
for i, ds  in enumerate(datasources):
    AA[i] = getPPInetwork(species, ds)

#homology correction
homInd = np.where(datasources == 'homology')[0]
tmInd = np.where(datasources == 'textmining')[0]
coInd = np.where(datasources == 'cooccurence')[0]

AA[tmInd] = AA[tmInd] * (1 - AA[homInd])
AA[coInd] = AA[coInd] * (1 - AA[homInd])

'''
#remove datasources that won't pass the threshold
tobedel = []
for i, (A, ds) in enumerate(zip(AA, datasources)):
    if np.max(A) <= 0. or (thresholdType == 'fixed' and np.max(A) < threshold):
        tobedel.append(i)

AA = np.delete(AA, tobedel, axis=0)
datasources = np.delete(datasources, tobedel)
'''

#number of possible combinations
F = len(datasources)
ss = 0
for i in range(1, F + 1):
	ss += binom(F, i)

ss = int(ss)
n_folds = 5

cv = KFold(n_splits=n_folds, shuffle=True, random_state=656391)

np.random.seed(1901273)

directory = '../experiments/' + species + '_' + annoPrefix + thresholdType

try:
    os.mkdir(directory)
except FileExistsError:
    pass

for i in range(n_folds):
    try:
        os.mkdir(directory + '/fold' + str(i))
    except FileExistsError:
        pass

try:
    os.mkdir(directory + '/networks')
except FileExistsError:
    pass


termNames = np.array(termNames)
for fold, (train, test) in enumerate(cv.split(Y)):

    Ytrain = Y[train]
    Ytest = Y[test]

    with open(directory + '/fold' + str(fold) + '/train.names', 'w') as fw:
        for g in geneNames[train]:
            fw.write(g + '\n')

    with open(directory + '/fold' + str(fold) + '/test.names', 'w') as fw:
        for g in geneNames[test]:
            fw.write(g + '\n')

    nonempty = np.where(np.sum(Ytest, 0) > 0)[0]
    Ytest = Ytest[:, nonempty]
    Ytrain = Ytrain[:, nonempty]

    termNames2 = termNames[nonempty]

    nonempty = np.where(np.sum(Ytrain, 0) > 0)[0]
    Ytest = Ytest[:, nonempty]
    Ytrain = Ytrain[:, nonempty]
    termNames2 = list(termNames2[nonempty])

    assert np.min(np.sum(Y, 1)) > 0

    with open(directory + '/fold' + str(fold) + '/terms.names', 'w') as fw:
        for g in termNames2:
            fw.write(g + '\n')

    #parentsCoord = getParentsCoord(termNames2, 'P', dag, mapping)
    #ic = calculateIC(Ytrain, parentsCoord)

    #np.save(directory + '/fold' + str(fold) + '/icvec.npy', ic)


try:
    os.mkdir('../data/' + species + '/networks/')
except FileExistsError:
    pass


counter = 0
for i in range(F + 1):

	if i == 0:
		A = Aexp

		ii, jj = np.where(A)
		with open('../data/' + species + '/networks/tmp' + str(counter) + '.network', 'w') as fw:
			for ind1, ind2 in zip(ii, jj):
				if ind1 > ind2:
					fw.write(str(ind1) + '\t' + str(ind2) + '\n')

	else:

		for j,p in enumerate(combinations(range(F), i)):
			counter += 1
			print(counter)
			sys.stdout.flush()
			Apred = integrateStringScores(AA[list(p)])

			if thresholdType == 'median':
				threshold = np.median(Apred[Apred > 0])


			Apred = (Apred >= threshold).astype(int)


			A = np.maximum(Aexp, Apred)

			ii, jj = np.where(A)
			with open('../data/' + species + '/networks/tmp' + str(counter) + '.network', 'w') as fw:
				for ind1, ind2 in zip(ii, jj):
					if ind1 > ind2:
						fw.write(str(ind1) + '\t' + str(ind2) + '\n')


	#print('starting n2v')
	#cmd = 'python node2vec/src/main.py --input tmp.network --output ../networks/emb' + str(i) + '.emb'
	#call(cmd, shell=True)
	#break
