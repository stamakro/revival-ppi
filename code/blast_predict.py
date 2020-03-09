from scipy.sparse import csr_matrix, load_npz
import pickle
import sys
import numpy as np
from sklearn.metrics import roc_auc_score
from routinesgo import *

dag, mapping = read_ontology_from_file('../data/go/go-final.obo')

experimentPath = sys.argv[1]
stuff = experimentPath.split('/')[2]
stuff2 = stuff.split('_')
species = stuff2[0]

if 'noIPI' in stuff:
    annoPrefix = stuff2[1] + '_' + stuff2[2] + '_'
else:
    annoPrefix = stuff2[1] + '_'


with open('../blast/annotations/' + annoPrefix + 'geneNames.pkl', 'rb') as f:
    geneNames = pickle.load(f)

gene2row = dict()
for i, g in enumerate(geneNames):
    gene2row[g] = i


with open('../blast/annotations/' + annoPrefix + 'termNames.pkl', 'rb') as f:
    termNames = pickle.load(f)

term2col = dict()
for i, t in enumerate(termNames):
    term2col[t] = i

Y = load_npz('../blast/annotations/' + annoPrefix + 'Y.npz')
parentsCoord = getParentsCoord(list(termNames), 'P', dag, mapping)
ic = calculateIC(Y.toarray(), parentsCoord)
freq = np.sum(Y.toarray(),0) / Y.shape[0]

icDict = dict()
freqDict = dict()
for t, ict, freqt in zip(termNames, ic, freq):
    icDict[t] = ict
    freqDict[t] = freqt

with open(experimentPath + 'icDict.pkl', 'wb') as f:
    pickle.dump(icDict, f)

with open(experimentPath + 'freqDict.pkl', 'wb') as f:
    pickle.dump(freqDict, f)

nfolds = 5

for fold in range(nfolds):

    testNames = set()
    testNamesList = []
    with open(experimentPath + 'fold' + str(fold) + '/test.names') as f:
        for line in f:
            testNames.add(line[:-1])
            testNamesList.append(line[:-1])

    predictions = dict()

    with open('../blast/' + species + '_fold' + str(fold) + '.txt') as f:
        for line in f:
            if line[0] != '#':
                fields = line.split()

                query = fields[0].split('|')[1]
                hit = fields[1].split('|')[1]
                identity = float(fields[2]) / 100.

                assert query in testNames

                evalue = float(fields[-2])

                if evalue < 0.01:
                    ypredTemp = Y[gene2row[hit]].toarray() * identity

                    if query in predictions:
                        predictions[query] = np.maximum(ypredTemp, predictions[query])
                    else:
                        predictions[query] = ypredTemp


    coverage = len(predictions) / len(testNames)
    Ytrue = np.zeros((len(testNames), Y.shape[1]))
    Ypred = np.zeros((len(testNames), Y.shape[1]))


    for i, g in enumerate(testNamesList):
        Ytrue[i] = Y[gene2row[g]].toarray()
        try:
            Ypred[i] = predictions[g]
        except KeyError:
            Ypred[i] = np.zeros((Y.shape[1]))


    tobedelTerms = np.where(np.sum(np.maximum(Ytrue, Ypred), 0) == 0)[0]

    Ytrue = np.delete(Ytrue, tobedelTerms, axis=1)
    Ypred = np.delete(Ypred, tobedelTerms, axis=1)

    newTermNames = np.delete(termNames, tobedelTerms)
    newIC = np.delete(ic, tobedelTerms)
    newFreq = np.delete(freq, tobedelTerms)

    #np.save(experimentPath + 'fold' + str(fold) + '/termFreq.npy', newFreq)
    #np.save(experimentPath + 'fold' + str(fold) + '/icvec.npy', newIC)


    posteriorDict = dict()

    posteriorDict['terms'] = newTermNames

    for i, g in enumerate(testNamesList):
        posteriorDict[g] = {'ytrue': Ytrue[i], 'ypred': Ypred[i]}

    with open(experimentPath + 'fold' + str(fold) + '/blast.pkl', 'wb') as f:
        pickle.dump(posteriorDict, f)
