from sklearn.model_selection import KFold
import numpy as np
from scipy.stats import pearsonr, spearmanr
import sys
from routinesnetworks import *

species = sys.argv[1]


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


n_folds = 5

degree = np.zeros((Aexp.shape[0]))
fmax = np.zeros((Aexp.shape[0]))
smin = np.zeros((Aexp.shape[0]))

cv = KFold(n_splits=n_folds, shuffle=True, random_state=656391)

np.random.seed(1901273)

current = 0
for fold, (train, test) in enumerate(cv.split(Aexp)):

    N = test.shape[0]

    degree[current:current+N] = np.sum(Aexp[test][:, train], 1)



    with open('../experiments/' + species + '_P_median/fold' + str(fold) + '/performance_gba_1hop.pkl', 'rb') as f:
        _, fmax_fold, _,  smin_fold, _ = pickle.load(f)

    fmax[current:current+N] = fmax_fold
    smin[current:current+N] = smin_fold


    current += N


#predictions = {'terms': termNames, 'gt': csr_matrix(Ytest)}
