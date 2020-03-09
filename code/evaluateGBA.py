import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score
import sys
import os
from routinesnetworks import evaluate

experimentPath = sys.argv[1]
nhops = int(sys.argv[2])
currentFold = int(sys.argv[3])

assert nhops == 1 or nhops == 2
nFolds = 5

classifiers = ['blast', 'biogrid', 'blast+biogrid']

for _, _, files in os.walk(experimentPath + 'fold0/ppi'):
    nString = len(files)

total_fmax = np.zeros((len(classifiers) + nString - 1,))
total_smin = np.zeros((len(classifiers) + nString - 1,))
total_nsmin = np.zeros((len(classifiers) + nString - 1,))


local_fmax = np.zeros((len(classifiers) + nString - 1,))
local_smin = np.zeros((len(classifiers) + nString - 1,))
local_nsmin = np.zeros((len(classifiers) + nString - 1,))

coverage = np.zeros((len(classifiers) + nString - 1,))

troc = np.zeros((len(classifiers) + nString - 1,))
pauc = np.zeros((len(classifiers) + nString - 1,))

with open(experimentPath + 'icDict.pkl', 'rb') as f:
    icDict = pickle.load(f)

for fold in range(nFolds):
    if fold != currentFold:
        continue

    print(fold)
    fullPath = experimentPath + 'fold' + str(fold) + '/'
    testNames = set()
    with open(fullPath + 'test.names') as f:
        for line in f:
            testNames.add(line[:-1])

    with open(fullPath + 'blast.pkl', 'rb') as f:
        blastResults = pickle.load(f)

    termNames = blastResults['terms']
    del blastResults['terms']
    ic = np.array([icDict[tt] for tt in termNames])
    term2col = dict()

    for jj, tt in enumerate(termNames):
        term2col[tt] = jj

    Ytrue = np.zeros((len(blastResults), len(termNames)))
    Ypost = np.zeros((len(blastResults), len(termNames)))
    for i, p in enumerate(blastResults):
        Ytrue[i] = blastResults[p]['ytrue']
        Ypost[i] = blastResults[p]['ypred']

    assert np.min(np.sum(Ytrue, 1)) > 0
    pred = np.where(np.max(Ypost, axis=1) > 0.)[0]

    coverage[0] = pred.shape[0] / len(testNames)

    total_fmax[0], total_smin[0], total_nsmin[0] = evaluate(Ytrue, Ypost, ic)

    Ytrue_short = Ytrue[pred]
    Ypost_short = Ypost[pred]

    local_fmax[0], local_smin[0], local_nsmin[0] = evaluate(Ytrue_short, Ypost_short, ic)

    #--------------------------------------------------------------------------
    with open(fullPath + 'ppi/0_knn.pkl', 'rb') as f:
        biogrid = pickle.load(f)

    ic2 = np.array([icDict[tt] for tt in biogrid['terms']])

    term2colBiogrid = dict()
    for kk, tt in enumerate(biogrid['terms']):
        term2colBiogrid[tt] = kk

    Ytrue = biogrid['gt'].toarray()
    Ypost_1h = biogrid['0'][nhops-1].toarray()

    assert np.min(np.sum(Ytrue, 1)) > 0
    pred = np.unique(np.where(np.logical_not(np.isnan(Ypost_1h)))[0])

    coverage[1] = pred.shape[0] / len(testNames)

    total_fmax[1], total_smin[1], total_nsmin[1] = evaluate(Ytrue, Ypost_1h, ic2)

    Ytrue2 = Ytrue[pred]
    Ypost2_1h = Ypost_1h[pred]
    local_fmax[1], local_smin[1], local_nsmin[1] = evaluate(Ytrue2, Ypost2_1h, ic2)

    print('combo\'s', flush=True)

    allterms = list(set(termNames).union(set(biogrid['terms'])))
    ic3 = np.zeros((len(allterms),))
    term2colFull = dict()
    for kk, tt in enumerate(allterms):
        term2colFull[tt] = kk
        ic3[kk] = icDict[tt]

    Ytrue_combo = np.zeros((len(testNames), len(allterms)))
    Ypred_combo1 = np.ones((len(testNames), len(allterms)))

    Ypost_1h[np.isnan(Ypost_1h)] = 0.

    for jj, tt in enumerate(allterms):
        if tt in term2colBiogrid:
            Ytrue_combo[:, jj] = Ytrue[:, term2colBiogrid[tt]]
            Ypred_combo1[:, jj] *= (1 - Ypost_1h[:, term2colBiogrid[tt]])


        if tt in term2col:
            Ypred_combo1[:, jj] *= (1 - Ypost[:, term2col[tt]])


    Ypred_combo1 = 1 - Ypred_combo1


    pred = np.where(np.max(Ypred_combo1, axis=1) > 0)[0]
    coverage[2] = pred.shape[0] / len(testNames)
    #print('Coverage Biogrid 2-hop: %3.1f' % (100 * pc_coverage), end='%\n')

    total_fmax[2], total_smin[2], total_nsmin[2] = evaluate(Ytrue_combo, Ypred_combo1, ic3)

    Ytrue2 = Ytrue_combo[pred]
    Ypost2_1h = Ypred_combo1[pred]

    local_fmax[2], local_smin[2], local_nsmin[2] = evaluate(Ytrue2, Ypost2_1h, ic3)

    for stringDs in range(1, 512):
        print('fold %d, network %d' % (fold, stringDs), flush=True)


        with open(fullPath + 'ppi/' + str(stringDs) + '_knn.pkl', 'rb') as f:
            biogrid = pickle.load(f)

        #Ytrue = biogrid['gt'].toarray()
        Ypost_1h = biogrid[nhops-1].toarray()

        assert np.min(np.sum(Ytrue, 1)) > 0
        pred = np.unique(np.where(np.logical_not(np.isnan(Ypost_1h)))[0])

        coverage[len(classifiers) + stringDs - 1] = pred.shape[0] / len(testNames)
        #print('Coverage Biogrid 1-hop: %3.1f' % (100 * pc_coverage), end='%\n')

        total_fmax[len(classifiers) + stringDs - 1], total_smin[len(classifiers) + stringDs - 1], total_nsmin[len(classifiers) + stringDs - 1] = evaluate(Ytrue, Ypost_1h, ic2)

        Ytrue2 = Ytrue[pred]
        Ypost2_1h = Ypost_1h[pred]
        local_fmax[len(classifiers) + stringDs - 1], local_smin[len(classifiers) + stringDs - 1], local_nsmin[len(classifiers) + stringDs - 1] = evaluate(Ytrue2, Ypost2_1h, ic2)

with open(fullPath + 'performance_gba_' + str(nhops) + 'hop.pkl', 'wb') as f:
    pickle.dump((coverage, total_fmax, local_fmax, total_smin, local_smin), f)
