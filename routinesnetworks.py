import numpy as np
from sklearn.metrics import precision_recall_fscore_support, average_precision_score, roc_auc_score, f1_score
import pickle
from scipy.stats import rankdata, spearmanr, pearsonr
from scipy.linalg import expm
import matplotlib.pyplot as plt
from sklearn.model_selection import KFold
import sys
from copy import deepcopy
from routinesgo import normalizedSemanticDistance

def predict(Atest, Ytrain, Ytest, ic, thresholds = None):

	if thresholds is None:
		thresholds = np.linspace(0.0, 1.0, 21)

	Y_posteriors = gbaPredict(Atest, Ytrain)

	pc_auprc = average_precision_score(Ytest.T, Y_posteriors.T, average=None)

	pc_f1 = np.zeros((thresholds.shape[0], Ytest.shape[0]))
	pc_nsd = np.zeros((thresholds.shape[0], Ytest.shape[0]))

	for i, thres in enumerate(thresholds):
		Ypred = (Y_posteriors >= thres).astype(int)

		pc_nsd[i] = normalizedSemanticDistance(Ytest, Ypred, ic)[2]
		pc_f1[i] = f1_score(Ytest.T, Ypred.T, average=None)


	return pc_auprc, pc_f1[np.argmax(np.mean(pc_f1, 1))], pc_nsd[np.argmin(np.mean(pc_nsd, 1))]


def predictBayes(Atest, Ytrain, Ytest, ic, thresholds = None):

	if thresholds is None:
		thresholds = np.linspace(0.0, 1.0, 21)


	priors = np.sum(Ytrain, 0) / float(Ytrain.shape[0])
	assert priors.shape[0] == Ytrain.shape[1]


	Y_likelihood = gbaPredict(Atest, Ytrain)


	marginal = Y_likelihood * priors + (1 - Y_likelihood) * (1 - priors)

	Y_posteriors = Y_likelihood * priors / marginal

	pc_auprc = average_precision_score(Ytest.T, Y_posteriors.T, average=None)

	pc_f1 = np.zeros((thresholds.shape[0], Ytest.shape[0]))
	pc_nsd = np.zeros((thresholds.shape[0], Ytest.shape[0]))

	for i, thres in enumerate(thresholds):
		Ypred = (Y_posteriors >= thres).astype(int)

		pc_nsd[i] = normalizedSemanticDistance(Ytest, Ypred, ic)[2]
		pc_f1[i] = f1_score(Ytest.T, Ypred.T, average=None)


	return pc_auprc, pc_f1[np.argmax(np.mean(pc_f1, 1))], pc_nsd[np.argmin(np.mean(pc_nsd, 1))]

def gbaPredict(A, Ytrain):

	assert A.shape[1] == Ytrain.shape[0]

	#Ypred = np.zeros((A.shape[0], Ytrain.shape[1]))


	Ypred = (A.dot(Ytrain.astype(float)).T  / np.sum(A,1)).T

	assert Ypred.shape == (A.shape[0], Ytrain.shape[1])

	'''

	for i, a in enumerate(A):
		neighbors = np.where(a)[0]

		if len(neighbors) == 0:
			Ypred[i] = np.zeros((Ytrain.shape[1],))
		else:
			Ypred[i] = np.sum(Ytrain[neighbors], axis=0) / float(neighbors.shape[0])
	'''
	return Ypred

def knnPredict(A, Ytrain, k):

	assert A.shape[1] == Ytrain.shape[0]

	Ypred = np.zeros((A.shape[0], Ytrain.shape[1]))

	knns = np.argsort(A, axis=1)[:, -k:]

	Ypred =  np.sum(Ytrain[knns], axis=1) / float(k)

	return Ypred


def calculateTOM(A):

	Atemp = A.astype(float)
	L = Atemp.dot(Atemp)

	degree = np.sum(A, 0)

	kk = np.tile(degree, (degree.shape[0],1))

	K = np.minimum(kk, kk.T)
	return (L+A) / (K - A + 1.)



def matchNetworkAndLabels(A, ppiGene2row, ppiRow2gene, geneNames):

	#if we have interactions, but not functions --> remove from ppi adjacency matrix
	tobedel = np.sort([ppiGene2row[g] for g in ppiGene2row if g not in geneNames])
	geneOrderPPI = np.array([ppiRow2gene[i] for i in range(len(ppiRow2gene))])

	A = np.delete(A, tobedel, axis=0)
	A = np.delete(A, tobedel, axis=1)
	geneOrderPPI = np.delete(geneOrderPPI, tobedel)

	ppiGene2row = dict()
	ppiRow2gene = dict()
	for i, g in enumerate(geneOrderPPI):
		ppiGene2row[g] = i
		ppiRow2gene[i] = g

	#if we have functions, but not interactions --> add empty rows/columns to adjacency matrix
	noInteractions = [g for g in geneNames if g not in ppiGene2row]

	A = np.hstack((A,np.zeros((A.shape[0], len(noInteractions)), int)))
	A = np.vstack((A,np.zeros((len(noInteractions), A.shape[1]), int)))

	current = len(ppiGene2row)
	for g in noInteractions:

		ppiGene2row[g] = current
		ppiRow2gene[current] = g

		current += 1


	ind = np.zeros((len(geneNames),), int)

	for i,g in enumerate(geneNames):
		ind[i] = ppiGene2row[g]

	A = A[ind]
	A = A[:, ind]

	assert np.max(np.abs(A - A.T)) == 0
	assert A.shape[0] == len(geneNames)


	return A, ppiGene2row, ppiRow2gene


def goLoader(species):

	with open('../data/' + species + '/annotations/Y.pkl') as f:
		Y = pickle.load(f).toarray()

	with open('../data/' + species + '/annotations/geneNames.pkl') as f:
		geneNames = pickle.load(f)

	with open('../data/' + species + '/annotations/termNames.pkl') as f:
		termNames = pickle.load(f)

	gene2row = dict()
	for i, g in enumerate(geneNames):
		gene2row[g] = i

	term2col = dict()
	for i, t in enumerate(termNames):
		term2col[t] = i

	return [Y, geneNames, termNames, gene2row, term2col]



def getPPInetwork(species, datasource):

	if datasource == 'biogrid':

		with open('../data/' + species + '/interactions/biogrid-final/A.pkl') as f:
			A = pickle.load(f).toarray()

		with open('../data/' + species + '/interactions/biogrid-final/row2protein.pkl') as f:
			ppiRow2gene = pickle.load(f)


	elif datasource == 'dl':
		print('gene & row info not available, but equal to result of matchNetworkAndLabels')
		with open('../data/' + species + '/interactions/dl/A.pkl') as f:
			A = pickle.load(f).toarray()

		A += A.T

		A = (A >= 0.5).astype(float)

		return [A, None, None]

	else:
		with open('../data/' + species + '/interactions/string-final/' + datasource + '.pkl') as f:
			A = pickle.load(f).toarray()

		with open('../data/' + species + '/interactions/biogrid-final/row2protein.pkl') as f:
			ppiRow2gene = pickle.load(f)

	A += A.T

	ppiGene2row = dict()
	for k in ppiRow2gene:
		ppiGene2row[ppiRow2gene[k]] = k


	return [A, ppiGene2row, ppiRow2gene]


def getMultiplePPInetworks(species, datasources):

	if type(datasources) != np.ndarray:
		datasources = np.array(datasources)

	#see string db faq for details
	for i, ds in enumerate(datasources):
		with open('../data/' + species + '/interactions/string-final/' + ds + '.pkl') as f:
			At = pickle.load(f).toarray()

		if i == 0:
			A = np.zeros((datasources.shape[0], At.shape[0], At.shape[1]))

		A[i] = At + At.T


	#homology correction
	if 'homology' in datasources:
		homInd = np.where(datasources == 'homology')[0]


		if 'textmining' in datasources:
			tmInd = np.where(datasources == 'textmining')[0]

			A[tmInd] = A[tmInd] * (1 - A[homInd])

		if 'cooccurence' in datasources:
			coInd = np.where(datasources == 'cooccurence')[0]
			A[coInd] = A[coInd] * (1 - A[homInd])


	Atotal = integrateStringScores(A)

	with open('../data/' + species + '/interactions/biogrid-final/row2protein.pkl') as f:
		ppiRow2gene = pickle.load(f)

	ppiGene2row = dict()
	for k in ppiRow2gene:
		ppiGene2row[ppiRow2gene[k]] = k


	return [Atotal, ppiGene2row, ppiRow2gene]



def integrateStringScores(As, p = 0.041):
	if np.max(As) > 1.:
		As /= 1000.

	#trick to ignore 0's
	As[As == 0.0] = p

	nopr = (As - p) / (1 - p)

	running = 1 - np.prod(1 - nopr, axis=0)

	final = running  + p *(1 - running)

	#reverse trick
	final[(final - p)**2 < 1e-6] = 0.0

	return final



def TestIntegrateStringScores(scores, p = 0.041):
	print('Do not use me')
	scores_noprior = (scores - p) / (1. - p)

	runningTotal = 1. - np.prod(1 - scores_noprior)

	scoresTotal = runningTotal + p * (1 - runningTotal)

	return scoresTotal


def compareTwoNets(Aexp, Apred, threshold = None):

	if threshold is None:
		threshold = np.median(Apred[Apred > 0])

	Apred = (Apred >= threshold).astype(int)

	Acomb = np.maximum(Apred, Aexp)

	diff = Acomb - Aexp

	stats = np.zeros((3,))
	stats[0] = np.sum(diff) / 2.0
	#print('Number of new interactions: %d' % (np.sum(diff) / 2.0))

	degreeExp = np.sum(Aexp, 0)
	degreeComb = np.sum(Acomb, 0)

	previouslyEmpty = np.where(np.logical_and(degreeExp == 0, degreeComb > 0))[0]

	#print('Proteins that had no partners, but do now: %d' % previouslyEmpty.shape[0])
	stats[1] = previouslyEmpty.shape[0]

	(iis, jjs) = np.where(diff)

	degreeRanks = rankdata(degreeExp)


	edgeHubness = np.zeros((np.sum(diff) / 2, ), float)
	c = -1
	for i,j in zip(iis, jjs):
		if i > j:
			c += 1
			#print(i,j)
			if degreeExp[i] == 0 or degreeExp[j] == 0:
				stats[2] += 1

	return stats




def exportToCytoscape(filename, Acomb, Aexp = None):
	if filename[-4:] != '.sif':
		filename += '.sif'
	with open(filename, 'w') as f:
		for i in xrange(Acomb.shape[0]-1):
			for j in xrange(i+1, Acomb.shape[0]):
				if Acomb[i,j]:

					if Aexp is not None:
						if Aexp[i,j]:
							f.write(str(i) + ' exp ' + str(j) + '\n')
						else:
							f.write(str(i) + ' pred ' + str(j) + '\n')

					else:
						f.write(str(i) + ' interaction ' + str(j) + '\n')



def cvLoop(A, Y, nfolds, classifier='gba', seed=656391, **kwargs):

	resultsDictionary = dict()

	if classifier == 'tom':
		A = calculateTOM(A)
		k = kwargs['k']

	elif classifier	== 'gba':
		try:
			nhops = kwargs['nhops']
		except KeyError:
			nhops = 1

		if np.isinf(nhops):
			A = expm(A)

		elif nhops > 1:
			Atemp = deepcopy(A)
			for i in range(nhops-1):
				Atemp = Atemp.dot(A)

			A = Atemp
			del Atemp



	cv = KFold(n_splits=nfolds, shuffle=True, random_state=seed)

	perfs = np.zeros((2, nfolds))
	corrs = np.zeros((nfolds))

	termsWithoutPrediction = np.zeros((nfolds,))
	proteinsMissed = np.zeros((nfolds,))
	termsWorseThanRandom = np.zeros((nfolds,))

	for fold, (train, test) in enumerate(cv.split(A)):
		print (fold)
		sys.stdout.flush()
		Atest = A[test][:, train]

		Ytrain = Y[train]
		Ytest = Y[test]

		if classifier == 'tom':
			Y_posteriors = knnPredict(Atest, Ytrain,k)

		elif classifier == 'gba':
			Y_posteriors = gbaPredict(Atest, Ytrain)

		else:
			print ('Error! Unknown classifier')
			return None

		pc_auprc = average_precision_score(Ytest.T, Y_posteriors.T, average=None)

		nonempty = np.where(np.sum(Ytest, 0) > 0)[0]

		tc_auroc = roc_auc_score(Ytest[:, nonempty], Y_posteriors[:, nonempty], average=None)

		perfs[0, fold] = np.nanmean(pc_auprc)
		perfs[1, fold] = np.nanmean(tc_auroc)

		tmpP = pc_auprc[np.logical_not(np.isnan(pc_auprc))]
		tmpD = np.sum(Atest, 1)[np.logical_not(np.isnan(pc_auprc))]

		corrs[fold] = spearmanr(tmpP, tmpD)[0]

		termsWithoutPrediction[fold] = np.sum(np.max(Y_posteriors, 0)[nonempty] == 0) / float(nonempty.shape[0])

		proteinsMissed[fold] = np.sum(np.isnan(pc_auprc)) / float(Atest.shape[0])
		termsWorseThanRandom[fold] = np.sum(tc_auroc <= 0.5) / float(nonempty.shape[0])

	resultsDictionary['missedProteins'] = proteinsMissed
	resultsDictionary['missedTerms'] = termsWithoutPrediction
	resultsDictionary['badTerms'] = termsWorseThanRandom


	resultsDictionary['auprc_p'] = perfs[0]
	resultsDictionary['auroc_t'] = perfs[1]
	resultsDictionary['perf-neighbors'] = corrs

	#print(np.mean(perfs), np.std(perfs))

	return resultsDictionary



def removeEdges(A, frac, new=True):
	#frac: the fraction of missing edges


	if frac == 1.0:
		return np.zeros(A.shape, int)

	if frac == 0.0:
		return A

	edges = np.where(A)


	edgePairs = [(i,j) for i, j in zip(edges[0], edges[1]) if i < j]

	toremove = int(np.round(len(edgePairs) * frac))

	pairsToRemove = np.random.permutation(edgePairs)[:toremove]

	if not new:

		for p in pairsToRemove:
			A[p[0], p[1]] = 0
			A[p[1], p[0]] = 0


		return A

	Anew = deepcopy(A)
	for p in pairsToRemove:
		Anew[p[0], p[1]] = 0
		Anew[p[1], p[0]] = 0


	return Anew
