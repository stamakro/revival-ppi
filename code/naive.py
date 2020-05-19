import pickle
from routinesnetworks import goLoader, evaluate
import numpy as np
import sys
from sklearn.model_selection import KFold

experimentPath = sys.argv[1]
stuff = experimentPath.split('/')[2]
stuff2 = stuff.split('_')
species = stuff2[0]

if 'noIPI' in stuff:
    annoPrefix = stuff2[1] + '_' + stuff2[2] + '_'
else:
    annoPrefix = stuff2[1] + '_'

thresholdType = stuff2[2]

print('Loading go...')
[Y, geneNames, termNames, gene2row, term2col] = goLoader(species)

with open(experimentPath + 'freqDict.pkl', 'rb') as f:
    freqDict = pickle.load(f)

with open(experimentPath + 'icDict.pkl', 'rb') as f:
    icDict = pickle.load(f)

freq = np.array([freqDict[t] for t in termNames])
ic = np.array([icDict[t] for t in termNames])

n_folds = 5
cv = KFold(n_splits=n_folds, shuffle=True, random_state=656391)
np.random.seed(1901273)

termNames = np.array(termNames)

fmax = np.zeros((n_folds,))
smin = np.zeros((n_folds,))

for fold, (train, test) in enumerate(cv.split(Y)):
    print(fold)

    Ytest = Y[test]
    Ypred = np.tile(freq, (Ytest.shape[0], 1))

    fmax[fold], smin[fold], _ = evaluate(Ytest, Ypred, ic, np.linspace(0, 1.,101))

print(np.mean(fmax), np.std(fmax, ddof=1))


with open(experimentPath + 'naive.pkl', 'wb') as f:
    pickle.dump([(np.mean(fmax), np.std(fmax, ddof=1)), (np.mean(smin), np.std(smin, ddof=1))], f)
