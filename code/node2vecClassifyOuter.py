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
from scipy.stats import ttest_rel
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


nwalksMatrix = [5, 10, 15, 20, 30]
pMatrix = [0.1, 0.5, 1.0, 1.5, 3.0]
qMatrix = [0.1, 0.5, 1.0, 1.5, 3.0]
lengthMatrix = [20, 40, 60, 80, 100, 160]
dimsMatrix = [32, 64, 96, 128, 150, 200, 500]
kMatrix = [1, 2, 5, 7, 9 , 11, 15, 21, 31, 51]
lMatrix = [0.0001, 0.001, 0.01, 0.1, 0.5, 1., 2., 3., 5., 10.]
assert len(kMatrix) == len(lMatrix)
clfMatrix = ['knn', 'lr']

performance = np.zeros((n_folds, len(nwalksMatrix), len(pMatrix), len(qMatrix), len(lengthMatrix), len(dimsMatrix), len(kMatrix), len(clfMatrix)))

for i1 in range(n_folds):
    for i2 in range(len(nwalksMatrix)):
        for i3 in range(len(pMatrix)):
            for i4 in range(len(qMatrix)):
                for i5 in range(len(lengthMatrix)):
                    with open('../experiments/n2v/' + species + '/' + str(i1)+ '_' + str(nwalksMatrix[i2]) + '_' + str(pMatrix[i3]) + '_' + str(qMatrix[i4]) + '_' + str(lengthMatrix[i5]) + '.txt') as f:
                        for line in f:
                            info, knn, lr = line.split('\t')
                            i67 = info.split(',')
                            i6 = int(i67[0][4:])
                            i7 = int(i67[1][11:])

                            performance[i1, i2, i3, i4, i5, i6, i7, 0] = float(knn[4:])
                            performance[i1, i2, i3, i4, i5, i6, i7, 1] = float(lr[3:])

nwalks = []
p = []
q = []
walkLength = []
dims = []
ks = []
clfs = []

for fold in range(n_folds):
    bestInd = np.argmax(performance[fold])
    nwi, pi, qi, wli, di, ki, ci = np.unravel_index(bestInd, performance[fold].shape)

    nwalks.append(nwalksMatrix[nwi])
    p.append(pMatrix[pi])
    q.append(qMatrix[qi])

    walkLength.append(lengthMatrix[wli])
    dims.append(dimsMatrix[di])

    if ci == 0:
        clfs.append('knn')
        ks.append(kMatrix[ki])

    else:
        clfs.append('lr')
        ks.append(lMatrix[ki])

print(clfs)
print(ks)
#sys.exit(0)

fmax = np.zeros((n_folds,3))
precision = np.zeros((n_folds,3))
recall = np.zeros((n_folds,3))
ru = np.zeros((n_folds,3))
mi = np.zeros((n_folds,3))
smin = np.zeros((n_folds,3))
coverage = np.zeros((n_folds, 3))



for fold, (train, test) in enumerate(cv.split(Y)):
    print(fold)


    Ytrain = Y[train]
    Ytest = Y[test]


    fileName = 'node2vec/emb/' + species + '_' + str(nwalks[fold]) + '_' + str(p[fold]) + '_' + str(q[fold]) + '_' + str(walkLength[fold]) + '_' + str(dims[fold]) + '.emb'
    print(fileName)
    #continue

    with open(fileName) as f:
        for i, line in enumerate(f):
            if i == 0:
                n2vDim = int(line.split()[1])
                assert n2vDim == dims[fold]
                X = np.zeros((Y.shape[0], n2vDim))
            else:
                fields = line.split()
                assert len(fields) == n2vDim + 1
                X[int(fields[0])] = np.array(fields[1:]).astype(float)


    Xtrain = X[train]
    Xtest = X[test]

    iits = np.where(np.sum(np.abs(Xtrain), 1) > 0)[0]

    Xtrain2 = Xtrain[iits]
    Ytrain2 = Ytrain[iits]



    if clfs[fold] == 'knn':
        k = ks[fold]
        nnk = np.argsort(cdist(Xtest, Xtrain2, metric='cosine'), axis=1)[:, :k]


        ypred = np.sum(Ytrain2[nnk], axis=1) / k

    else:
        w = ridgeTrain(Xtrain2, Ytrain2, ks[fold])
        ypred = ridgePredict(Xtest, w)


    empty = np.where(np.sum(np.abs(Xtest), 1) == 0)[0]
    ypred[empty] = 0.0

    coverage[fold, 0] = np.mean(np.max(ypred,1)>0)

    with open(experimentPath +  'fold' + str(fold) + '/node2vec/biogrid.pkl', 'wb') as f:
        pickle.dump(ypred, f)



    #fmax[fold, 0], smin[fold,0], _ = evaluate(Ytest, ypred, ic)
    precision[fold, 0], fmax[fold, 0], ru[fold, 0], mi[fold,0], smin[fold,0] = evaluateFull(Ytest, ypred, ic)

'''
performance = np.zeros((n_folds, len(nwalksMatrix), len(pMatrix), len(qMatrix), len(lengthMatrix), len(dimsMatrix), len(kMatrix), len(clfMatrix)))

for i1 in range(n_folds):
    for i2 in range(len(nwalksMatrix)):
        for i3 in range(len(pMatrix)):
            for i4 in range(len(qMatrix)):
                for i5 in range(len(lengthMatrix)):
                    with open('../experiments/n2v/' + species + '/string/' + str(i1)+ '_' + str(nwalksMatrix[i2]) + '_' + str(pMatrix[i3]) + '_' + str(qMatrix[i4]) + '_' + str(lengthMatrix[i5]) + '.txt') as f:
                        for line in f:
                            info, knn, lr = line.split('\t')
                            i67 = info.split(',')
                            i6 = int(i67[0][4:])
                            i7 = int(i67[1][11:])

                            performance[i1, i2, i3, i4, i5, i6, i7, 0] = float(knn[4:])
                            performance[i1, i2, i3, i4, i5, i6, i7, 1] = float(lr[3:])

nwalks = []
p = []
q = []
walkLength = []
dims = []
ks = []
clfs = []

for fold in range(n_folds):
    bestInd = np.argmax(performance[fold])
    nwi, pi, qi, wli, di, ki, ci = np.unravel_index(bestInd, performance[fold].shape)

    nwalks.append(nwalksMatrix[nwi])
    p.append(pMatrix[pi])
    q.append(qMatrix[qi])

    walkLength.append(lengthMatrix[wli])
    dims.append(dimsMatrix[di])

    if ci == 0:
        clfs.append('knn')
        ks.append(kMatrix[ki])

    else:
        clfs.append('lr')
        ks.append(lMatrix[ki])


if species == 'yeast':
    stringNumber = 39

elif species == 'ecoli':
    stringNumber = 129

elif species == 'arabidopsis':
    stringNumber = 39
elif species == 'tomato':
    stringNumber = 88


print(clfs)
print(ks)





for fold, (train, test) in enumerate(cv.split(Y)):
    print(fold)


    Ytrain = Y[train]
    Ytest = Y[test]


    fileName = 'node2vec/emb/string/' + species + '_' + str(nwalks[fold]) + '_' + str(p[fold]) + '_' + str(q[fold]) + '_' + str(walkLength[fold]) + '_' + str(dims[fold]) + '.emb'

    with open(fileName) as f:
        for i, line in enumerate(f):
            if i == 0:
                n2vDim = int(line.split()[1])
                assert n2vDim == dims[fold]
                X = np.zeros((Y.shape[0], n2vDim))
            else:
                fields = line.split()
                assert len(fields) == n2vDim + 1
                X[int(fields[0])] = np.array(fields[1:]).astype(float)


    Xtrain = X[train]
    Xtest = X[test]

    iits = np.where(np.sum(np.abs(Xtrain), 1) > 0)[0]

    Xtrain2 = Xtrain[iits]
    Ytrain2 = Ytrain[iits]


    if clfs[fold] == 'knn':
        k = ks[fold]
        nnk = np.argsort(cdist(Xtest, Xtrain2, metric='cosine'), axis=1)[:, :k]


        ypred = np.sum(Ytrain2[nnk], axis=1) / k

    else:
        w = ridgeTrain(Xtrain2, Ytrain2, ks[fold])
        ypred = ridgePredict(Xtest, w)


    empty = np.where(np.sum(np.abs(Xtest), 1) == 0)[0]
    ypred[empty] = 0.0

    coverage[fold, 1] = np.mean(np.max(ypred,1)>0)


    fmax[fold,1], smin[fold,1], _ = evaluate(Ytest, ypred, ic)

    with open(experimentPath +  'fold' + str(fold) + '/node2vec/biogrid_string.pkl', 'wb') as f:
        pickle.dump(ypred, f)
'''

performance = np.zeros((n_folds, len(nwalksMatrix), len(pMatrix), len(qMatrix), len(lengthMatrix), len(dimsMatrix), len(kMatrix), len(clfMatrix)))

for i1 in range(n_folds):
    for i2 in range(len(nwalksMatrix)):
        for i3 in range(len(pMatrix)):
            for i4 in range(len(qMatrix)):
                for i5 in range(len(lengthMatrix)):
                    with open('../experiments/n2v/' + species + '/withUnannotated/' + str(i1)+ '_' + str(nwalksMatrix[i2]) + '_' + str(pMatrix[i3]) + '_' + str(qMatrix[i4]) + '_' + str(lengthMatrix[i5]) + '.txt') as f:
                        for line in f:
                            info, knn, lr = line.split('\t')
                            i67 = info.split(',')
                            i6 = int(i67[0][4:])
                            i7 = int(i67[1][11:])

                            performance[i1, i2, i3, i4, i5, i6, i7, 0] = float(knn[4:])
                            performance[i1, i2, i3, i4, i5, i6, i7, 1] = float(lr[3:])

nwalks = []
p = []
q = []
walkLength = []
dims = []
ks = []
clfs = []

for fold in range(n_folds):
    bestInd = np.argmax(performance[fold])
    nwi, pi, qi, wli, di, ki, ci = np.unravel_index(bestInd, performance[fold].shape)

    nwalks.append(nwalksMatrix[nwi])
    p.append(pMatrix[pi])
    q.append(qMatrix[qi])

    walkLength.append(lengthMatrix[wli])
    dims.append(dimsMatrix[di])

    if ci == 0:
        clfs.append('knn')
        ks.append(kMatrix[ki])

    else:
        clfs.append('lr')
        ks.append(lMatrix[ki])



print(clfs)
print(ks)





for fold, (train, test) in enumerate(cv.split(Y)):
    print(fold)


    Ytrain = Y[train]
    Ytest = Y[test]


    fileName = 'node2vec/emb/withUnannotated/' + species + '_' + str(nwalks[fold]) + '_' + str(p[fold]) + '_' + str(q[fold]) + '_' + str(walkLength[fold]) + '_' + str(dims[fold]) + '.emb'
    print(fileName)
    #continue
    with open(fileName) as f:
        for i, line in enumerate(f):
            if i == 0:
                n2vDim = int(line.split()[1])
                assert n2vDim == dims[fold]
                X = np.zeros((Y.shape[0], n2vDim))
            else:
                fields = line.split()
                assert len(fields) == n2vDim + 1
                try:
                    X[int(fields[0])] = np.array(fields[1:]).astype(float)
                except IndexError:
                    pass

    Xtrain = X[train]
    Xtest = X[test]

    iits = np.where(np.sum(np.abs(Xtrain), 1) > 0)[0]

    Xtrain2 = Xtrain[iits]
    Ytrain2 = Ytrain[iits]


    if clfs[fold] == 'knn':
        k = ks[fold]
        nnk = np.argsort(cdist(Xtest, Xtrain2, metric='cosine'), axis=1)[:, :k]


        ypred = np.sum(Ytrain2[nnk], axis=1) / k

    else:
        w = ridgeTrain(Xtrain2, Ytrain2, ks[fold])
        ypred = ridgePredict(Xtest, w)


    empty = np.where(np.sum(np.abs(Xtest), 1) == 0)[0]
    ypred[empty] = 0.0

    coverage[fold, 2] = np.mean(np.max(ypred,1)>0)


    #fmax[fold,2], smin[fold,2], _ = evaluate(Ytest, ypred, ic)
    precision[fold, 2], fmax[fold, 2], ru[fold, 2], mi[fold,2], smin[fold,2] = evaluateFull(Ytest, ypred, ic)
print('Fmax')
print(np.mean(fmax,0), np.std(fmax,0, ddof=1))
'''print('\nPrecision')
print(np.mean(precision,0), np.std(precision,0, ddof=1))'''
print('\n\nSmin')
print(np.mean(smin,0), np.std(smin,0, ddof=1))
'''print('\nRU')
print(np.mean(ru,0), np.std(ru,0, ddof=1))
print('\nMI')
print(np.mean(mi,0), np.std(mi,0, ddof=1))'''
print('\n\nCoverage')
print(np.mean(coverage,0), np.std(coverage,0, ddof=1))

print(species)
print('without unannotated - with unannotated\n\n')

tt = ttest_rel(fmax[:,0], fmax[:, 2])
print('Fmax: %f\t%f' % (tt[0], tt[1]))

tt = ttest_rel(smin[:,0], smin[:, 2])
print('Smin: %f\t%f' % (tt[0], tt[1]))

tt = ttest_rel(coverage[:,0], coverage[:, 2])
print('Cvrg: %f\t%f' % (tt[0], tt[1]))


