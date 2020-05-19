import numpy as np
import pickle
import sys
from sklearn.model_selection import KFold
from scipy.spatial.distance import cdist
from scipy.sparse import csr_matrix
from routinesnetworks import *
from routinesgo import semanticDistance
import os
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import f1_score
from copy import deepcopy
from scipy.special import expit
from datetime import datetime

def inverse(X):
    m, n = X.shape
    assert m == n
    return np.linalg.solve(X, np.eye(m))


def ridgeTrain(X, Y, lamda):
    Y[Y==0] = -1
    return inverse(X.T.dot(X) + lamda * np.eye(X.shape[1])).dot(X.T).dot(Y)

def ridgePredict(X, W):
    return expit(X.dot(W))


experimentPath = sys.argv[1]
stuff = experimentPath.split('/')[2]
stuff2 = stuff.split('_')
species = stuff2[0]

if 'noIPI' in stuff:
    annoPrefix = stuff2[1] + '_' + stuff2[2] + '_'
else:
    annoPrefix = stuff2[1] + '_'

thresholdType = stuff2[2]


[Y, geneNames, termNames, gene2row, term2col] = goLoader(species)

with open(experimentPath + 'icDict.pkl', 'rb') as f:
    icDict = pickle.load(f)

ic = np.array([icDict[t] for t in termNames])


n_folds = 5


cv = KFold(n_splits=n_folds, shuffle=True, random_state=656391)

np.random.seed(1901273)

termNames = np.array(termNames)

currentFold = int(sys.argv[2])
nwalks = int(sys.argv[3])
p = float(sys.argv[4])
q = float(sys.argv[5])
walkLength = int(sys.argv[6])


try:
    isitString = int(sys.argv[7])
except IndexError:
    isitString = False



if isitString:
    outfile = '../experiments/n2v/' + species + '/string/' + str(currentFold) + '_' + str(nwalks) + '_' + str(p) + '_' + str(q) + '_' + str(walkLength) + '.txt'
else:
    outfile = '../experiments/n2v/' + species + '/' + str(currentFold) + '_' + str(nwalks) + '_' + str(p) + '_' + str(q) + '_' + str(walkLength) + '.txt'


dims = [32, 64, 96, 128, 150, 200, 500]
ks = [1, 2, 5, 7, 9, 11, 15, 21, 31, 51]
lamdas = [0.0001, 0.001, 0.01, 0.1, 0.5, 1., 2., 3., 5., 10.]


fmaxKnn = np.zeros((len(dims), len(ks)))
fmaxLR = np.zeros((len(dims), len(lamdas)))


with open(outfile, 'w') as fw:

    for fold, (train, test) in enumerate(cv.split(Y)):

        if fold != currentFold:
            continue

        Ntrain = int(train.shape[0] * 0.8)
        permTrain = np.random.permutation(train)

        train_in = permTrain[:Ntrain]
        valid_in = permTrain[Ntrain:]

        Ytrain = Y[train_in]
        Yval = Y[valid_in]


        for idim, dim in enumerate(dims):
            if isitString:
                fileName = 'node2vec/emb/string/' + species + '_' + str(nwalks) + '_' + str(p) + '_' + str(q) + '_' + str(walkLength) + '_' + str(dim) + '.emb'
            else:
                fileName = 'node2vec/emb/' + species + '_' + str(nwalks) + '_' + str(p) + '_' + str(q) + '_' + str(walkLength) + '_' + str(dim) + '.emb'

            with open(fileName) as f:
                for i, line in enumerate(f):
                    if i == 0:
                        n2vDim = int(line.split()[1])
                        assert n2vDim == dim
                        X = np.zeros((Y.shape[0], n2vDim))
                    else:
                        fields = line.split()
                        assert len(fields) == n2vDim + 1
                        X[int(fields[0])] = np.array(fields[1:]).astype(float)


            Xtrain = X[train_in]
            Xval = X[valid_in]

            iitr = np.where(np.sum(np.abs(Xtrain), 1) > 0)[0]
            iits = np.where(np.sum(np.abs(Xval), 1) > 0)[0]

            Xtrain = Xtrain[iitr]
            Ytrain_temp = Ytrain[iitr]
            Xval = Xval[iits]
            Yval_temp = Yval[iits]

            Ypred = np.zeros(Yval.shape)

            nn = np.argsort(cdist(Xval, Xtrain, metric='cosine'), axis=1)#[:, :5]

            for i, (k, l) in enumerate(zip(ks, lamdas)):
                print('%d, %d' % (idim, i), flush=True)

                nnk = nn[:,:k]
                ypred = np.sum(Ytrain_temp[nnk], axis=1) / k

                fmaxKnn[idim, i], _, _ = evaluate(Yval_temp, ypred, ic)

                w = ridgeTrain(Xtrain, Ytrain_temp, l)

                fmaxLR[idim, i], _, _ = evaluate(Yval_temp, ridgePredict(Xval, w), ic)

                fw.write('Dim ' + str(idim) + ', parameter ' + str(i) + '\tknn ' + str(fmaxKnn[idim, i]) + '\tlr ' + str(fmaxLR[idim, i]) + '\n')
