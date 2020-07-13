import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pickle
import sys
from scipy.special import binom
from itertools import combinations

def index2ds(ind, datasources):
	if ind == 0:
		return 'exp'

	Nds = datasources.shape[0]
	if ind > 2 ** Nds - 1:

		raise ValueError

	counter = 1
	for i in range(1,ds):
		iterations = int(binom(ds, i))
		if counter + iterations > ind:

			for j, combo in enumerate(combinations(range(Nds), i)):
				if counter + j == ind:
					combo = list(combo)
					return datasources[combo]

		counter += iterations

	return 'all'


def ds2ind(datasources):
	#d2i = dict()
	Nds = datasources.shape[0]

	d2i = np.zeros((Nds, 2**Nds), bool)

	counter = 1
	for i in range(1,ds):
		iterations = int(binom(ds, i))
		for combo in combinations(range(Nds), i):
			combo = list(combo)
			d2i[combo, counter] = 1

			counter += 1

	d2i[:, -1] = 1
	return d2i



species = sys.argv[1]
experimentPath = '../experiments/' + species + '_P_median/'

if species == 'tomato':
	datasources = np.array(['neighborhood_transferred', 'coexpression_transferred', 'experiments_transferred', 'textmining', 'cooccurence','fusion', 'textmining_transferred', 'homology'])

elif species == 'ecoli':
	datasources = np.array(['neighborhood', 'neighborhood_transferred', 'coexpression_transferred', 'experiments_transferred', 'textmining', 'cooccurence','fusion', 'textmining_transferred', 'homology'])


else:
	datasources = np.array(['coexpression', 'neighborhood_transferred', 'coexpression_transferred', 'experiments_transferred', 'textmining', 'cooccurence','fusion', 'textmining_transferred', 'homology'])


metrics = ['coverage', 'global_fmax', 'local_fmax', 'global_smin', 'local_smin']

ds = datasources.shape[0]

d2i = ds2ind(datasources)

ind_hom = np.where(d2i[np.where(datasources == 'homology')[0]])[1]

ind_tm = np.where(np.logical_or(d2i[np.where(datasources == 'textmining')[0]], d2i[np.where(datasources == 'textmining_transferred')[0]]))[1]

with open(experimentPath + 'naive.pkl', 'rb') as f:
	naive = pickle.load(f)

clf = 'gba_1hop'
n_folds = 5
for fold in range(n_folds):
    with open(experimentPath + 'fold' + str(fold) + '/performance_' + clf + '.pkl', 'rb') as f:
        data = pickle.load(f)

    if fold == 0:
        bigMatrix = np.zeros((n_folds, len(metrics), data[0].shape[0]))

    for i, v in enumerate(data):
        bigMatrix[fold, i] = v

meanMatrix = np.mean(bigMatrix, axis=0)
stdMatrix = np.std(bigMatrix, axis=0, ddof=1)

meanMatrix = np.delete(meanMatrix, 0, axis=0)
stdMatrix = np.delete(stdMatrix, 0, axis=0)

blast = meanMatrix[:,0]
biogrid = meanMatrix[:,1]
blastAndBiogrid = meanMatrix[:,2]

meanMatrix = np.delete(meanMatrix, [0, 2], axis=1)
stdMatrix = np.delete(stdMatrix, [0,2], axis=1)

stringStart = 1


mS = np.zeros((ds+1,))
mF = np.zeros((ds+1,))
fig = plt.figure(figsize=(12,9))

ax_g_f = fig.add_subplot(1, 2, 1)
ax_g_s = fig.add_subplot(1, 2, 2)

counter = stringStart

xxf = np.zeros(meanMatrix.shape[1])

for i in range(ds+1):
    if i == 0:
        ax_g_f.scatter(0, biogrid[0], color='C0', edgecolor='k')
        ax_g_f.axvline(1, color='k', linestyle=':', alpha=0.3)

        ax_g_s.scatter(0, biogrid[2], color='C0', edgecolor='k')
        ax_g_s.axvline(1, color='k', linestyle=':', alpha=0.3)
        mF[0] = biogrid[0]
        mS[0] = biogrid[2]


    else:
        iterations = int(binom(ds, i))

        mF[i] = np.mean(meanMatrix[0, counter:counter+iterations])
        mS[i] = np.mean(meanMatrix[2, counter:counter+iterations])
        xxf[counter:counter+iterations] = 2*i + 1.5 * np.random.rand(iterations) - 0.75

        ax_g_f.scatter(xxf[counter:counter+iterations], meanMatrix[0, counter:counter+iterations], color='C0', edgecolor='k' )
        ax_g_f.axvline(2*i+1, color='k', linestyle=':', alpha=0.3)

        ax_g_s.scatter(xxf[counter:counter+iterations], meanMatrix[2, counter:counter+iterations], color='C0', edgecolor='k' )
        ax_g_s.axvline(2*i+1, color='k', linestyle=':', alpha=0.3)


        counter += iterations



ax_g_f.scatter(xxf[ind_tm], meanMatrix[0, ind_tm], color='y', edgecolor='k', label='text mining')
ax_g_f.scatter(xxf[ind_hom], meanMatrix[0, ind_hom], color='k', edgecolor='k', label='homology')
ax_g_f.scatter(xxf[np.intersect1d(ind_tm,ind_hom)], meanMatrix[0, np.intersect1d(ind_tm,ind_hom)], color='k', edgecolor='y', label='text mining + homology')
ax_g_f.scatter([], [], color='C0', edgecolor='k', label='rest')

ax_g_s.scatter(xxf[ind_tm], meanMatrix[2, ind_tm], color='y', edgecolor='k', label='text mining')
ax_g_s.scatter(xxf[ind_hom], meanMatrix[2, ind_hom], color='k', edgecolor='k', label='homology')
ax_g_s.scatter(xxf[np.intersect1d(ind_tm,ind_hom)], meanMatrix[2, np.intersect1d(ind_tm,ind_hom)], color='k', edgecolor='y', label='text mining + homology')
ax_g_s.scatter([], [], color='C0', edgecolor='k', label='rest')

ax_g_f.plot(np.arange(0,2*mF.shape[0], 2), mF, color='C1')
ax_g_s.plot(np.arange(0, 2*mS.shape[0], 2), mS, color='C1')
ax_g_f.set_xticks([2*i for i in range(ds+2)])
ax_g_f.set_xticklabels([str(i) for i in range(ds+2)])
ax_g_s.set_xticks([2*i for i in range(ds+2)])
ax_g_s.set_xticklabels([str(i) for i in range(ds+2)])
ax_g_f.set_xlabel('# of data sources', fontsize=20)
ax_g_f.set_ylabel('Fmax', fontsize=20)
ax_g_s.set_xlabel('# of data sources', fontsize=20)
ax_g_s.set_ylabel('Smin', fontsize=20)
ax_g_f.set_xlim(-1, 2*ds+2)
ax_g_s.set_xlim(-1, 2*ds+2)

if clf == 'gba_1hop' and species == 'yeast':
    ax_g_f.set_ylim(0.3, 0.6)

ax_g_f.axhline(naive[0][0], color='k', label='naive')
ax_g_s.axhline(naive[1][0], color='k', label='naive')
ax_g_f.axhline(blast[1], color='C3', label='blast')
ax_g_s.axhline(blast[3], color='C3', label='blast')
ax_g_f.axhline(blastAndBiogrid[0], color='C2', linestyle='--', label='blast+biogrid')
ax_g_s.axhline(blastAndBiogrid[2], color='C2', linestyle='--', label='blast+biogrid')

ax_g_f.legend(prop={'size':12})
ax_g_s.legend(prop={'size':12})
fig.tight_layout()



fig.savefig('../figures/main-paper/fig4_' + species + '.png')
fig.savefig('../figures/main-paper/fig4_' + species + '.eps')


