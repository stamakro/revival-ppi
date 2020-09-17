import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pickle
import sys
from routinesnetworks import evaluate
from scipy.stats import ttest_rel

experimentPath = sys.argv[1]
stuff = experimentPath.split('/')[2]
stuff2 = stuff.split('_')
species = stuff2[0]

fignames = ['coverage', 'global', 'local']

clf = 'gba_1hop'
n_folds = 5

if species == 'tomato':
	ds = 8
else:
	ds = 9

metrics = ['coverage', 'global_fmax', 'local_fmax', 'global_smin', 'local_smin']

Ncombos = 2**ds - 1
with open(experimentPath + 'icDict.pkl', 'rb') as f:
    icDict = pickle.load(f)

for fold in range(n_folds):
    with open(experimentPath + 'fold' + str(fold) + '/performance_' + clf + '.pkl', 'rb') as f:
        data = pickle.load(f)

    if fold == 0:
        bigMatrix = np.zeros((n_folds, len(metrics), data[0].shape[0]))
        bigMatrixWeighted = np.zeros((n_folds, len(metrics), Ncombos))

    for i, v in enumerate(data):
        bigMatrix[fold, i] = v

    fullPath = experimentPath + 'fold' + str(fold) + '/'
    with open(fullPath + 'ppi/0_knn.pkl', 'rb') as f:
        biogrid = pickle.load(f)
        Ytrue = biogrid['gt'].toarray()
        terms = biogrid['terms']

    ic2 = np.array([icDict[tt] for tt in biogrid['terms']])


    for i in range(1, Ncombos+1):
        print('Fold %d, combo %d/%d' % (fold, i, Ncombos), flush=True)
        with open(experimentPath + 'fold' + str(fold) + '/ppi/' + str(i) +  '_knn_weightedString.pkl', 'rb') as f:
            Ypost = pickle.load(f)[0].toarray()

        pred = np.unique(np.where(np.logical_not(np.isnan(Ypost)))[0])

        bigMatrixWeighted[fold, 0, i-1] = pred.shape[0] / Ytrue.shape[0]

        Ypost[np.where(np.isnan(Ypost))] = 0.0

        bigMatrixWeighted[fold, 1, i-1], bigMatrixWeighted[fold, 2, i-1], _ = evaluate(Ytrue, Ypost, ic2)

    assert np.max(bigMatrixWeighted[fold]) > 0


with open('../results/weightedString_' + species + '.pkl', 'wb') as f:
	pickle.dump({'normal': bigMatrix, 'weighted': bigMatrixWeighted}, f)

metrics = ['Coverage', 'Fmax', 'Smin']
'''

with open('../results/weightedString_' + species + '.pkl', 'rb') as f:
	results = pickle.load(f)

bigMatrix = results['normal']
bigMatrixWeighted = results['weighted']

'''
bigMatrix = np.delete(bigMatrix, [0, 1, 2], axis=2)
bigMatrix = np.delete(bigMatrix, [2, 4], axis=1)


bigMatrixWeighted = np.delete(bigMatrixWeighted, [3, 4], axis=1)


meanMatrix = np.mean(bigMatrix, axis=0)
stdMatrix = np.std(bigMatrix, axis=0, ddof=1)


meanMatrixWeighted = np.mean(bigMatrixWeighted, axis=0)
stdMatrixWeighted = np.std(bigMatrixWeighted, axis=0, ddof=1)


iis = [1, 2, 0]

fig = plt.figure()

for i, ii in enumerate(iis):

	ax = fig.add_subplot(1, 3, i+1, aspect='auto')
	ax.scatter(meanMatrix[ii], meanMatrixWeighted[ii], color='C0', edgecolor='k')

	m = np.minimum(np.min(meanMatrix[ii]), np.min(meanMatrixWeighted[ii]))
	M = np.maximum(np.max(meanMatrix[ii]), np.max(meanMatrixWeighted[ii]))

	xx = np.linspace(0.85*m, 1.15*M)
	ax.plot(xx, xx, 'k--')
	ax.set_title(metrics[ii])

	ax.set_xlabel('Binary, 50% highest edges', fontsize=9)
	ax.set_ylabel('Weighted, all edges', fontsize=9)

	print(metrics[ii])
	
	if ii == 2:
		optB = np.argmin(meanMatrix[ii])	
		optW = np.argmin(meanMatrixWeighted[ii])	
	else:
		optB = np.argmax(meanMatrix[ii])	
		optW = np.argmax(meanMatrixWeighted[ii])	

	print('binary: %.2f +/- %.3f\nweighted: %.2f +/- %.3f' % (meanMatrix[ii,optB], stdMatrix[ii, optB], meanMatrixWeighted[ii, optW], stdMatrixWeighted[ii, optW]))
	print(ttest_rel(bigMatrix[:, ii, optB], bigMatrixWeighted[:, ii, optW]))
	print('\n\n')


plt.tight_layout()
fig.savefig('../figures/weightedString_' + species + '.png', dpi=600)
fig.savefig('../figures/weightedString_' + species + '.eps', dpi=1200)
