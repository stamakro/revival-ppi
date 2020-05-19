import numpy as np
import pickle
import sys
from itertools import combinations
from sklearn.model_selection import KFold
from scipy.special import binom
from scipy.spatial.distance import cdist
from scipy.sparse import csr_matrix
from routinesnetworks import *
from routinesgo import semanticDistance
import os
from sklearn.linear_model import LogisticRegression
from sklearn.multioutput import MultiOutputClassifier
from sklearn.metrics import f1_score


experimentPath = sys.argv[1]
stuff = experimentPath.split('/')[2]
stuff2 = stuff.split('_')
species = stuff2[0]

if 'noIPI' in stuff:
    annoPrefix = stuff2[1] + '_' + stuff2[2] + '_'
else:
    annoPrefix = stuff2[1] + '_'

thresholdType = stuff2[2]

classifier = sys.argv[2]
edgeFile = sys.argv[3]
nrEdgeFiles = int(sys.argv[4])


if len(sys.argv) == 5:
	missingFraction = 0.0
else:
	missingFraction = float(sys.argv[5])

print('Loading go...')
[Y, geneNames, termNames, gene2row, term2col] = goLoader(species)

with open(experimentPath + 'icDict.pkl', 'rb') as f:
    icDict = pickle.load(f)

ic = np.array([icDict[t] for t in termNames])

print('Loading experimental PPI network...')

Aexp = getPPInetwork(species, 'biogrid')
assert np.max(np.abs((Aexp - Aexp.T)) < 1e-10)

try:
    Aexp2 = getPPInetwork(species, 'experiments')
    assert np.max(np.abs((Aexp2 - Aexp2.T)) < 1e-10)

    if np.max(Aexp2) > 0.:
    	threshold = np.median(Aexp2[Aexp2 > 0])
    	Aexp2 = (Aexp2 >= threshold).astype(int)
    	Aexp = np.maximum(Aexp, Aexp2)
except FileNotFoundError:
    pass

Aexp = Aexp.astype(float)


degree = np.sum(Aexp,0)
with open('../data/' + species + '/degrees/0.pkl', 'wb') as f:
    pickle.dump(degree, f)



if species == 'tomato':
	datasources = np.array(['neighborhood_transferred', 'coexpression_transferred', 'experiments_transferred', 'textmining', 'cooccurence','fusion', 'textmining_transferred', 'homology'])

elif species == 'ecoli':
	datasources = np.array(['neighborhood', 'neighborhood_transferred', 'coexpression_transferred', 'experiments_transferred', 'textmining', 'cooccurence','fusion', 'textmining_transferred', 'homology'])


else:
	datasources = np.array(['coexpression', 'neighborhood_transferred', 'coexpression_transferred', 'experiments_transferred', 'textmining', 'cooccurence','fusion', 'textmining_transferred', 'homology'])


AA = np.zeros((datasources.shape[0], Aexp.shape[0], Aexp.shape[1]))
for i, ds in enumerate(datasources):
	print ('Reading predicted %d' % i)
	At = getPPInetwork(species, ds)
	AA[i] = At

print('Correcting for homology...')

#homology correction
homInd = np.where(datasources == 'homology')[0]

tmInd = np.where(datasources == 'textmining')[0]

AA[tmInd] = AA[tmInd] * (1 - AA[homInd])

coInd = np.where(datasources == 'cooccurence')[0]
AA[coInd] = AA[coInd] * (1 - AA[homInd])

Astring = integrateStringScores(AA[:])


if thresholdType == 'median':
	threshold = np.median(Astring[Astring > 0])

Astring = (Astring >= threshold).astype(float)



print('Loading predicted edges...')
Adl = np.zeros(Aexp.shape)
trind = np.triu_indices_from(Adl, k=1)
NN = trind[0].shape[0]
post = np.zeros(NN)

start = 0

for i in range(nrEdgeFiles):
    with open(edgeFile + str(i) + '.pkl', 'rb') as f:
        post_subset = pickle.load(f)
    N = post_subset.shape[0]
    post[start:start+N] = post_subset

    start += N

Adl[trind] = post
Adl += Adl.T

#threshold = np.median(post)
threshold = 0.5

Adl = (Adl > threshold).astype(float)
Adl[np.diag_indices_from(Adl)] = 0

Acombo = np.maximum(Aexp, Adl)
#Acombo = np.maximum(Acombo, Astring)

Agmd = np.maximum(Aexp, Astring)
ii = np.where(np.sum(Agmd,0) == 0)[0]
Agmd[ii] = Adl[ii]

n_folds = 5


cv = KFold(n_splits=n_folds, shuffle=True, random_state=656391)

np.random.seed(1901273)

termNames = np.array(termNames)

coverage = np.zeros((n_folds,2))
global_fmax = np.zeros((n_folds,2))
global_smin = np.zeros((n_folds,2))
local_fmax = np.zeros((n_folds,2))
local_smin = np.zeros((n_folds,2))

for fold, (train, test) in enumerate(cv.split(Y)):

    print(fold)

    Ytrain = Y[train]
    Ytest = Y[test]
    predictions = {'terms': termNames, 'gt': csr_matrix(Ytest)}

    if classifier == 'gba':
        pp_1hop = predict('gba', atest=Agmd[test][:, train], ytrain=Ytrain)
        pp_1hop_combo = predict('gba', atest=Acombo[test][:, train], ytrain=Ytrain)


        predind = np.where(np.max(pp_1hop, 1) > 0)[0]
        coverage[fold, 0] = predind.shape[0] / Ytest.shape[0]

        predind1 = np.where(np.max(pp_1hop_combo, 1) > 0)[0]
        coverage[fold, 1] = predind1.shape[0] / Ytest.shape[0]

        global_fmax[fold,0], global_smin[fold,0], _ = evaluate(Ytest, pp_1hop, ic)
        global_fmax[fold,1], global_smin[fold,1], _ = evaluate(Ytest, pp_1hop_combo, ic)

        if coverage[fold,0] < 1:

            local_fmax[fold,0], local_smin[fold,0], _ = evaluate(Ytest[predind], pp_1hop[predind], ic)

        else:
            local_fmax[fold,0], local_smin[fold,0] = global_fmax[fold,0], global_smin[fold,0]



        if coverage[fold,1] < 1:

            local_fmax[fold,1], local_smin[fold,1], _ = evaluate(Ytest[predind1], pp_1hop[predind1], ic)

        else:
            local_fmax[fold,1], local_smin[fold,1] = global_fmax[fold,1], global_smin[fold,1]


        #global_smin[counter] = semanticDistance(Ytest, Ypred, ic)[2]


    else:
        print('not implemented')


    '''
	elif classifier == 'n2v_lr' and not os.path.exists(experimentPath + 'fold' + str(fold_nr) + '/n2v/0_lr.pkl'):
		with open('../data/' + species + '/networks/tmp0.emb') as f:
			for i, line in enumerate(f):
				if i == 0:
					n2vDim = int(line.split()[1])
					X = np.zeros((Aexp.shape[0], n2vDim))
				else:
					fields = line.split()
					assert len(fields) == n2vDim + 1
					X[int(fields[0])] = np.array(fields[1:]).astype(float)
		Xtrain = X[train]
		Xtest = X[test]
		ii_tr = np.where(np.max(np.abs(Xtrain), axis=1) > 0.)[0]
		ii_ts = np.where(np.max(np.abs(Xtest), axis=1) > 0.)[0]
		clf = MultiOutputClassifier(LogisticRegression()).fit(Xtrain[ii_tr], Ytrain[ii_tr])
		y = clf.predict_proba(Xtest[ii_ts])

	elif classifier == 'n2v_knn' and not os.path.exists(experimentPath + 'fold' + str(fold_nr) + '/n2v/0_knn.pkl'):
		with open('../data/' + species + '/networks/tmp0.emb') as f:
			for i, line in enumerate(f):
				if i == 0:
					n2vDim = int(line.split()[1])
					X = np.zeros((Aexp.shape[0], n2vDim))
				else:
					fields = line.split()
					assert len(fields) == n2vDim + 1
					X[int(fields[0])] = np.array(fields[1:]).astype(float)
		Xtrain = X[train]
		Xtest = X[test]
		Ypred = np.zeros(Ytest.shape)
		ii_tr = np.where(np.max(np.abs(Xtrain), axis=1) > 0.)[0]
		ii_ts = np.where(np.max(np.abs(Xtest), axis=1) > 0.)[0]
		coverage[counter] = ii_ts.shape[0] / Xtest.shape[0]
		empty_ts = np.where(np.max(np.abs(Xtest), axis=1) == 0.)[0]
		k = 5
		nn = np.argsort(cdist(Xtest, Xtrain[ii_tr], metric='cosine'), axis=1)[:, :5]
		ypred = np.sum(Ytrain[nn], axis=1) / k
		ypred[empty_ts] = 0

		global_fmax[counter], global_smin[counter], _ = evaluate(Ytest, ypred, ic, np.linspace(0., 1., k+1))
		#global_smin[counter] = semanticDistance(Ytest, Ypred, ic)[2]

		if empty_ts.shape[0] > 0:
			local_fmax[counter], local_smin[counter], _ = evaluate(Ytest[ii_ts], ypred[ii_ts], ic, np.linspace(0., 1., k+1))

		else:
			local_fmax[counter] = global_fmax[counter]
			local_smin[counter] = global_smin[counter]
    '''
