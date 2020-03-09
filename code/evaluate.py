import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score
import sys
from routinesnetworks import evaluate

experimentPath = sys.argv[1]

nFolds = 5

classifiers = ['blast', 'biogrid_1hop', 'biogrid_2hop', 'blast+biogrid_1hop', 'blast+biogrid_2hop']

total_fmax = np.zeros((len(classifiers) + 511, nFolds))
total_smin = np.zeros((len(classifiers) + 511, nFolds))
total_nsmin = np.zeros((len(classifiers) + 511, nFolds))


local_fmax = np.zeros((len(classifiers) + 511, nFolds))
local_smin = np.zeros((len(classifiers) + 511, nFolds))
local_nsmin = np.zeros((len(classifiers) + 511, nFolds))

coverage = np.zeros((len(classifiers) + 511, nFolds))

troc = np.zeros((len(classifiers) + 511, nFolds))
pauc = np.zeros((len(classifiers) + 511, nFolds))

with open(experimentPath + 'icDict.pkl', 'rb') as f:
    icDict = pickle.load(f)

for fold in range(nFolds):
    print(fold)
    print(total_fmax)
    print('')
    
    print('')
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

    coverage[0, fold] = pred.shape[0] / len(testNames)
    #print('Coverage BLAST: %3.1f' % (100 * pc_coverage), end='%\n')

    total_fmax[0, fold], total_smin[0, fold], total_nsmin[0, fold] = evaluate(Ytrue, Ypost, ic)

    Ytrue_short = Ytrue[pred]
    Ypost_short = Ypost[pred]

    local_fmax[0, fold], local_smin[0, fold], local_nsmin[0, fold] = evaluate(Ytrue_short, Ypost_short, ic)

    #--------------------------------------------------------------------------
    with open(fullPath + 'ppi/0_knn.pkl', 'rb') as f:
        biogrid = pickle.load(f)

    ic2 = np.array([icDict[tt] for tt in biogrid['terms']])

    term2colBiogrid = dict()
    for kk, tt in enumerate(biogrid['terms']):
        term2colBiogrid[tt] = kk

    Ytrue = biogrid['gt'].toarray()
    Ypost_1h = biogrid['0'][0].toarray()

    assert np.min(np.sum(Ytrue, 1)) > 0
    pred = np.unique(np.where(np.logical_not(np.isnan(Ypost_1h)))[0])

    coverage[1, fold] = pred.shape[0] / len(testNames)
    #print('Coverage Biogrid 1-hop: %3.1f' % (100 * pc_coverage), end='%\n')

    total_fmax[1, fold], total_smin[1, fold], total_nsmin[1, fold] = evaluate(Ytrue, Ypost_1h, ic2)

    Ytrue2 = Ytrue[pred]
    Ypost2_1h = Ypost_1h[pred]
    local_fmax[1, fold], local_smin[1, fold], local_nsmin[1, fold] = evaluate(Ytrue2, Ypost2_1h, ic2)

    #---------------------------------------------------------------------------
    Ypost_2h = biogrid['0'][1].toarray()

    assert np.min(np.sum(Ytrue, 1)) > 0
    pred = np.unique(np.where(np.logical_not(np.isnan(Ypost_2h)))[0])

    coverage[2, fold] = pred.shape[0] / len(testNames)
    #print('Coverage Biogrid 2-hop: %3.1f' % (100 * pc_coverage), end='%\n')

    total_fmax[2, fold], total_smin[2, fold], total_nsmin[2, fold] = evaluate(Ytrue, Ypost_2h, ic2)

    Ytrue2 = Ytrue[pred]
    Ypost2_2h = Ypost_2h[pred]

    local_fmax[2, fold], local_smin[2, fold], local_nsmin[2, fold] = evaluate(Ytrue2, Ypost2_2h, ic2)

    print('combo\'s')

    allterms = list(set(termNames).union(set(biogrid['terms'])))
    ic3 = np.zeros((len(allterms),))
    term2colFull = dict()
    for kk, tt in enumerate(allterms):
        term2colFull[tt] = kk
        ic3[kk] = icDict[tt]

    Ytrue_combo = np.zeros((len(testNames), len(allterms)))
    Ypred_combo1 = np.ones((len(testNames), len(allterms)))
    Ypred_combo2 = np.ones((len(testNames), len(allterms)))

    Ypost_1h[np.isnan(Ypost_1h)] = 0.
    Ypost_2h[np.isnan(Ypost_2h)] = 0.

    for jj, tt in enumerate(allterms):
        if tt in term2colBiogrid:
            Ytrue_combo[:, jj] = Ytrue[:, term2colBiogrid[tt]]
            Ypred_combo1[:, jj] *= (1 - Ypost_1h[:, term2colBiogrid[tt]])
            Ypred_combo2[:, jj] *= (1 - Ypost_2h[:, term2colBiogrid[tt]])

        if tt in term2col:
            Ypred_combo1[:, jj] *= (1 - Ypost[:, term2col[tt]])
            Ypred_combo2[:, jj] *= (1 - Ypost[:, term2col[tt]])

    Ypred_combo1 = 1 - Ypred_combo1
    Ypred_combo2 = 1 - Ypred_combo2


    pred = np.where(np.max(Ypred_combo1, axis=1) > 0)[0]
    coverage[3, fold] = pred.shape[0] / len(testNames)
    #print('Coverage Biogrid 2-hop: %3.1f' % (100 * pc_coverage), end='%\n')

    total_fmax[3, fold], total_smin[3, fold], total_nsmin[3, fold] = evaluate(Ytrue_combo, Ypred_combo1, ic3)

    Ytrue2 = Ytrue_combo[pred]
    Ypost2_1h = Ypred_combo1[pred]

    local_fmax[3, fold], local_smin[3, fold], local_nsmin[3, fold] = evaluate(Ytrue2, Ypost2_1h, ic3)


    pred = np.where(np.max(Ypred_combo2, axis=1) > 0)[0]
    coverage[4, fold] = pred.shape[0] / len(testNames)
    #print('Coverage Biogrid 2-hop: %3.1f' % (100 * pc_coverage), end='%\n')

    total_fmax[4, fold], total_smin[4, fold], total_nsmin[4, fold] = evaluate(Ytrue_combo, Ypred_combo2, ic3)

    Ytrue2 = Ytrue_combo[pred]
    Ypost2_2h = Ypred_combo2[pred]

    local_fmax[4, fold], local_smin[4, fold], local_nsmin[4, fold] = evaluate(Ytrue2, Ypost2_2h, ic3)
    
    for stringDs in range(1, 512):
        print('fold %d, network %d' % (fold, stringDs))
        

        with open(fullPath + 'ppi/' + str(stringDs) + '_knn.pkl', 'rb') as f:
            biogrid = pickle.load(f)

        ic2 = np.array([icDict[tt] for tt in biogrid['terms']])

        term2colBiogrid = dict()
        for kk, tt in enumerate(biogrid['terms']):
            term2colBiogrid[tt] = kk

        Ytrue = biogrid['gt'].toarray()
        Ypost_1h = biogrid[str(stringDs)][0].toarray()

        assert np.min(np.sum(Ytrue, 1)) > 0
        pred = np.unique(np.where(np.logical_not(np.isnan(Ypost_1h)))[0])

        coverage[len(classifiers) + stringDs, fold] = pred.shape[0] / len(testNames)
        #print('Coverage Biogrid 1-hop: %3.1f' % (100 * pc_coverage), end='%\n')

        total_fmax[len(classifiers) + stringDs, fold], total_smin[len(classifiers) + stringDs, fold], total_nsmin[len(classifiers)+ stringDs, fold] = evaluate(Ytrue, Ypost_1h, ic2)

        Ytrue2 = Ytrue[pred]
        Ypost2_1h = Ypost_1h[pred]
        local_fmax[len(classifiers) + stringDs, fold], local_smin[len(classifiers) + stringDs, fold], local_nsmin[len(classifiers) + stringDs, fold] = evaluate(Ytrue2, Ypost2_1h, ic2)







