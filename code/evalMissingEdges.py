import numpy as np
import pickle
import sys
from sklearn.model_selection import KFold
from scipy.special import binom
from scipy.spatial.distance import cdist
from scipy.sparse import csr_matrix
from routinesnetworks import *
from routinesgo import semanticDistance
import os
from sklearn.metrics import f1_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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

fractions = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]

n_folds = 5
n_seeds = 5

local_fmax_R = np.zeros((n_seeds, len(fractions), n_folds))
total_fmax_R = np.zeros((n_seeds, len(fractions), n_folds))

local_smin_R = np.zeros((n_seeds, len(fractions), n_folds))
total_smin_R = np.zeros((n_seeds, len(fractions), n_folds))

coverage_R = np.zeros((n_seeds, len(fractions), n_folds))


local_fmax_D = np.zeros((n_seeds, len(fractions), n_folds))
total_fmax_D = np.zeros((n_seeds, len(fractions), n_folds))

local_smin_D = np.zeros((n_seeds, len(fractions), n_folds))
total_smin_D = np.zeros((n_seeds, len(fractions), n_folds))

coverage_D = np.zeros((n_seeds, len(fractions), n_folds))


termNames = np.array(termNames)


for seedIndex in range(n_seeds):

	for j, frac in enumerate(fractions):
		if j > 0:
			with open('../data/yeast/networks/downsampled/seed' + str(seedIndex) + '/random_' + str(frac) + '.pkl', 'rb') as f:
				Ar = pickle.load(f).toarray()

			with open('../data/yeast/networks/downsampled/seed' + str(seedIndex) + '/degree_' + str(frac) + '.pkl', 'rb') as f:
				Ad = pickle.load(f).toarray()

		cv = KFold(n_splits=n_folds, shuffle=True, random_state=656391)

		np.random.seed(1901273)
		for fold, (train, test) in enumerate(cv.split(Y)):
			print('%d %f %d' % (seedIndex, frac, fold))

			Ytrain = Y[train]
			Ytest = Y[test]
			predictions = {'terms': termNames, 'gt': csr_matrix(Ytest)}


			if j == 0:
				Ypost_1h = predict('gba', atest=Aexp[test][:, train], ytrain=Ytrain)
				
				pred = np.unique(np.where(np.logical_not(np.isnan(Ypost_1h)))[0])

				coverage_R[seedIndex, j, fold] = pred.shape[0] / Ytest.shape[0] 
				coverage_D[seedIndex, j, fold] = coverage_R[seedIndex, j, fold] 

				Ypost_1h[np.where(np.isnan(Ypost_1h))] = 0.0

				total_fmax_R[seedIndex, j, fold], total_smin_R[seedIndex, j, fold], _ = evaluate(Ytest, Ypost_1h, ic)
				total_fmax_D[seedIndex, j, fold] = total_fmax_R[seedIndex, j, fold]
				total_smin_D[seedIndex, j, fold] = total_smin_R[seedIndex, j, fold]

 
				Ytrue2 = Ytrain[pred]
				Ypost2_1h = Ypost_1h[pred]
				local_fmax_R[seedIndex, j, fold], local_smin_R[seedIndex, j, fold], _ = evaluate(Ytrue2, Ypost2_1h, ic)

				local_fmax_D[seedIndex, j, fold] = local_fmax_R[seedIndex, j, fold]
				local_smin_D[seedIndex, j, fold] = local_smin_R[seedIndex, j, fold]


			else:
				#random
				Ypost_1h = predict('gba', atest=Ar[test][:, train], ytrain=Ytrain)
				
				pred = np.unique(np.where(np.logical_not(np.isnan(Ypost_1h)))[0])

				coverage_R[seedIndex, j, fold] = pred.shape[0] / Ytest.shape[0]
 
				Ypost_1h[np.where(np.isnan(Ypost_1h))] = 0.0

				total_fmax_R[seedIndex, j, fold], total_smin_R[seedIndex, j, fold], _ = evaluate(Ytest, Ypost_1h, ic)

				Ytrue2 = Ytrain[pred]
				Ypost2_1h = Ypost_1h[pred]
				local_fmax_R[seedIndex, j, fold], local_smin_R[seedIndex, j, fold], _ = evaluate(Ytrue2, Ypost2_1h, ic)


				#degree
				Ypost_1h = predict('gba', atest=Ad[test][:, train], ytrain=Ytrain)
				
				pred = np.unique(np.where(np.logical_not(np.isnan(Ypost_1h)))[0])

				coverage_D[seedIndex, j, fold] = pred.shape[0] / Ytest.shape[0]
 
				Ypost_1h[np.where(np.isnan(Ypost_1h))] = 0.0

				total_fmax_D[seedIndex, j, fold], total_smin_D[seedIndex, j, fold], _ = evaluate(Ytest, Ypost_1h, ic)

				Ytrue2 = Ytrain[pred]
				Ypost2_1h = Ypost_1h[pred]
				local_fmax_D[seedIndex, j, fold], local_smin_D[seedIndex, j, fold], _ = evaluate(Ytrue2, Ypost2_1h, ic)



with open('../results/downsampling_yeast_random.pkl', 'wb') as f:
	pickle.dump({'c': coverage_R, 'gf': total_fmax_R, 'gs': total_smin_R, 'lf': local_fmax_R, 'ls': local_smin_R}, f)


with open('../results/downsampling_yeast_degree.pkl', 'wb') as f:
	pickle.dump({'c': coverage_D, 'gf': total_fmax_D, 'gs': total_smin_D, 'lf': local_fmax_D, 'ls': local_smin_D}, f)











