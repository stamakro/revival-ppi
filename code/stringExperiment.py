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
fold_nr = int(sys.argv[3])


if len(sys.argv) == 4:
	missingFraction = 0.0
else:
	missingFraction = float(sys.argv[4])


try:
	os.mkdir(experimentPath + 'fold' + str(fold_nr) + '/ppi/')
except FileExistsError:
	pass

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

#Aexp = removeEdges(Aexp, missingFraction)

F = len(datasources)
ss = 0
for i in range(1, F + 1):
	ss += binom(F, i)

ss = int(ss)

local_fmax = np.zeros((ss+1,))
global_fmax = np.zeros((ss+1,))

local_smin = np.zeros((ss+1,))
global_smin = np.zeros((ss+1,))
coverage = np.zeros((ss+1,))

n_folds = 5


cv = KFold(n_splits=n_folds, shuffle=True, random_state=656391)

np.random.seed(1901273)

termNames = np.array(termNames)

for fold, (train, test) in enumerate(cv.split(Y)):

	if fold != fold_nr:
		continue

	Ytrain = Y[train]
	Ytest = Y[test]
	predictions = {'terms': termNames, 'gt': csr_matrix(Ytest)}

	'''
	nonempty = np.where(np.sum(Ytest, 0) > 0)[0]
	Ytest = Ytest[:, nonempty]
	Ytrain = Ytrain[:, nonempty]

	termNames2 = termNames[nonempty]

	nonempty = np.where(np.sum(Ytrain, 0) > 0)[0]
	Ytest = Ytest[:, nonempty]
	Ytrain = Ytrain[:, nonempty]
	termNames2 = list(termNames2[nonempty])
	'''

	counter = 0
	for i in range(F + 1):


		if i == 0:
			if classifier == 'gba' and not os.path.exists(experimentPath + 'fold' + str(fold_nr) + '/ppi/0_knn.pkl'):
				pp_1hop = predict('gba', atest=Aexp[test][:, train], ytrain=Ytrain)
				
				A2 = Aexp.dot(Aexp)
				pp_2hop = predict('gba', atest=A2[test][:, train], ytrain=Ytrain)

				predictions['0'] = (csr_matrix(pp_1hop), csr_matrix(pp_2hop))

				with open(experimentPath + 'fold' + str(fold_nr) + '/ppi/0_knn.pkl', 'wb') as fw:
					pickle.dump(predictions, fw)

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


		else:

			for j,p in enumerate(combinations(range(F), i)):
				counter += 1

				#print('Fold:%d, Iteration %d/%d' % (fold, counter, ss))
				if classifier == 'gba' and os.path.exists(experimentPath + 'fold' + str(fold_nr) + '/ppi/' + str(counter) + '_knn.pkl'):
					continue
				print('Fold:%d, Iteration %d/%d' % (fold, counter, ss), flush=True)
				Apred = integrateStringScores(AA[list(p)])

				if thresholdType == 'median':
					threshold = np.median(Apred[Apred > 0])

				Apred = (Apred >= threshold).astype(float)

				if np.max(Apred) < 1:
					print('No new edges predicted for combination %d' % counter)
					continue

				A = np.maximum(Aexp, Apred)
				if classifier == 'gba':
					pp_1hop = predict('gba', atest=A[test][:, train], ytrain=Ytrain)

					A2 = A.dot(A)
					pp_2hop = predict('gba', atest=A2[test][:, train], ytrain=Ytrain)
					pp = (csr_matrix(pp_1hop), csr_matrix(pp_2hop))


					with open(experimentPath + 'fold' + str(fold_nr) + '/ppi/' + str(counter) + '_knn.pkl', 'wb') as fw:
						pickle.dump(pp, fw)

				elif classifier == 'n2v_knn':

					with open('../data/' + species + '/networks/tmp' + str(counter) + '.emb') as f:
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
					empty_ts = np.where(np.max(np.abs(Xtest), axis=1) == 0.)[0]
					coverage[counter] = ii_ts.shape[0] / Xtest.shape[0]
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

if classifier == 'n2v_knn':
	with open(experimentPath + 'fold' + str(fold_nr) + '/performance_n2v_5nn.pkl', 'wb') as f:
		pickle.dump((coverage, global_fmax, local_fmax, global_smin, local_smin), f)
