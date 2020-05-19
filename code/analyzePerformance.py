import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pickle
import sys
from scipy.special import binom
from scipy.stats import norm
from sklearn.mixture import GaussianMixture
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

def statisticalTest(performance, ind, metric, ptype='mean', Nperms=10000, seed=19912071):

	if metric == 'smin' and ptype == 'max':
		print('!!!warning!!!')

	fin = eval('np.' + ptype + '(performance[np.where(ind)[0]])')
	fout = eval('np.' + ptype + '(performance[np.where(ind==0)[0]])')

	if metric == 'fmax':
		diff = fin - fout
	else:
		diff = fout - fin

	rstats = np.zeros(Nperms)
	for i in range(Nperms):
		rind = np.random.permutation(ind)
		rfin = eval('np.' + ptype + '(performance[np.where(rind)[0]])')
		rfout = eval('np.' + ptype + '(performance[np.where(rind==0)[0]])')

		if metric == 'fmax':
			rstats[i] = rfin - rfout
		else:
			rstats[i] = rfout - rfin


	pvalue = np.mean(rstats >= diff)

	return fin, fout, diff, pvalue










experimentPath = sys.argv[1]
stuff = experimentPath.split('/')[2]
stuff2 = stuff.split('_')
species = stuff2[0]

#k = int(sys.argv[2])

clf = 'gba+blast_1hop'
clf = 'gba_1hop'
n_folds = 5

if species == 'tomato':
	datasources = np.array(['neighborhood_transferred', 'coexpression_transferred', 'experiments_transferred', 'textmining', 'cooccurence','fusion', 'textmining_transferred', 'homology'])

elif species == 'ecoli':
	datasources = np.array(['neighborhood', 'neighborhood_transferred', 'coexpression_transferred', 'experiments_transferred', 'textmining', 'cooccurence','fusion', 'textmining_transferred', 'homology'])


else:
	datasources = np.array(['coexpression', 'neighborhood_transferred', 'coexpression_transferred', 'experiments_transferred', 'textmining', 'cooccurence','fusion', 'textmining_transferred', 'homology'])


ds = datasources.shape[0]


metrics = ['coverage', 'global_fmax', 'local_fmax', 'global_smin', 'local_smin']

mS = np.zeros((ds+1,))
mF = np.zeros((ds+1,))
mS_local = np.zeros((ds+1,))
mF_local = np.zeros((ds+1,))
cov = np.zeros((ds+1,))


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

meanMatrix = np.delete(meanMatrix, [0, 2], axis=1)
stdMatrix = np.delete(stdMatrix, [0,2], axis=1)

biogrid = meanMatrix[:,0]
stringStart = 1

sys.exit(0)
xcoords = np.zeros(meanMatrix.shape[1])

counter = stringStart
for i in range(1,ds):
	iterations = int(binom(ds, i))
	xcoords[counter:counter+iterations] = 2*i + 1.5 * np.random.rand(iterations) - 0.75
	counter += iterations

xcoords[-1] = 2 * ds

y = meanMatrix[0].reshape(-1,1)
#clustering = KMeans(n_clusters=k).fit(meanMatrix[0].reshape(-1,1))

bics = []
aics = []
for kk in range(1, 7):
	gmm = GaussianMixture(n_components=kk, covariance_type='spherical', random_state=1991)
	gmm.fit(y)

	bics.append(gmm.bic(y))
	aics.append(gmm.aic(y))

kkBest = 1 + np.argmin(bics)

gmm = GaussianMixture(n_components=kkBest, covariance_type='spherical', random_state=1991)
gmm.fit(y)
labels = gmm.predict(y)


fig = plt.figure()
#0 coverage, 1 global, 2 local

colors=np.array(['C' + str(i) for i in range(kkBest)])
alphas = np.array([0.2, 1.])

ax_g_f = fig.add_subplot(1, 2, 1)
ax_g_s = fig.add_subplot(1, 2, 2)

ax_g_f.scatter(xcoords, meanMatrix[0], color=colors[labels], edgecolor='k')
ax_g_s.scatter(xcoords, meanMatrix[2], color=colors[labels], edgecolor='k')

for ll in range(1, 21, 2):
	ax_g_f.axvline(ll, color='k', linestyle=':', alpha=0.3)
	ax_g_s.axvline(ll, color='k', linestyle=':', alpha=0.3)

ax_g_f.set_xticks([2*i for i in range(ds+2)])
ax_g_f.set_xticklabels([str(i) for i in range(ds+2)])
ax_g_s.set_xticks([2*i for i in range(ds+2)])
ax_g_s.set_xticklabels([str(i) for i in range(ds+2)])
ax_g_f.set_xlabel('# of data sources', fontsize=20)
ax_g_f.set_ylabel('global Fmax', fontsize=20)
ax_g_s.set_xlabel('# of data sources', fontsize=20)
ax_g_s.set_ylabel('global Smin', fontsize=20)


fig.savefig('../figures/clustering/' + species + '_fg.png')


fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(y, bins=20, edgecolor='k', density=True)
xx = np.linspace(np.min(y)*0.9, np.max(y)*1.1, 200)

density = np.zeros(xx.shape)
for kk in range(kkBest):
	density += norm.pdf(xx, gmm.means_[kk], np.sqrt(gmm.covariances_[kk]))*gmm.weights_[kk]

ax.plot(xx, density, color='C1')
fig.savefig('../figures/clustering/' + species + '_fg_hist.png')

d2i = ds2ind(datasources)

pvals = np.zeros(datasources.shape[0])
fig = plt.figure(figsize=(16, 12))
print('Fmax')
print('%24s\t%s\t%s\t%s\t%s\t%s' % ('DATASOURCE', 'alone', 'in', 'out', 'diff', 'p-value'))
for i,d in enumerate(datasources):

	ax = fig.add_subplot(3,3, i+1)

	#ax.scatter(xcoords, meanMatrix[0], color=colors[d2i[i].astype(int)], s=12, alpha=alphas[d2i[i].astype(int)], edgecolor='k')
	ind1 = np.where(d2i[i])
	ind2 = np.setdiff1d(np.arange(meanMatrix.shape[1]), ind1)
	ax.scatter(xcoords[ind2], meanMatrix[0, ind2], color=colors[0], s=12, alpha=alphas[0], edgecolor='k', label='without')
	ax.scatter(xcoords[ind1], meanMatrix[0, ind1], color=colors[1], s=12, alpha=alphas[1], edgecolor='k', label='with')

	plt.legend(loc='lower right')
	fig.tight_layout()


	for ll in range(1, 21, 2):
		ax.axvline(ll, color='k', linestyle=':', alpha=0.3)

	ax.set_xticks([2*i for i in range(ds+2)])
	ax.set_xticklabels([str(i) for i in range(ds+2)])
	ax.set_xlabel('# of data sources')
	ax.set_ylabel('global Fmax')
	ax.set_title(d)


	falone = meanMatrix[0, 1 + i]
	fin, fout, diff, pvals[i] = statisticalTest(meanMatrix[0], d2i[i], 'fmax', ptype='mean')

	#pval = np.minimum(1., pvalue * datasources.shape[0])

	print('%24s\t%.3f\t%.3f\t%.3f\t%.4f\t%.5f' % (d, falone, fin, fout, diff, pvals[i]))
fig.savefig('../figures/clustering/per_data_source/' + species + '_fmax.png', dpi=500)

fig = plt.figure(figsize=(16, 12))
print('\n\nSmin')
print('%24s\t%s\t%s\t%s\t%s\t%s' % ('DATASOURCE', 'alone', 'in', 'out', 'diff', 'p-value'))
for i,d in enumerate(datasources):

	ax = fig.add_subplot(3,3, i+1)
	ind1 = np.where(d2i[i])
	ind2 = np.setdiff1d(np.arange(meanMatrix.shape[1]), ind1)
	ax.scatter(xcoords[ind2], meanMatrix[2, ind2], color=colors[0], s=12, alpha=alphas[0], edgecolor='k', label='without')
	ax.scatter(xcoords[ind1], meanMatrix[2, ind1], color=colors[1], s=12, alpha=alphas[1], edgecolor='k', label='with')
	plt.legend(loc='upper right')
	fig.tight_layout()

	for ll in range(1, 21, 2):
		ax.axvline(ll, color='k', linestyle=':', alpha=0.3)

	ax.set_xticks([2*i for i in range(ds+2)])
	ax.set_xticklabels([str(i) for i in range(ds+2)])
	ax.set_xlabel('# of data sources')
	ax.set_ylabel('global Smin')
	ax.set_title(d)

	salone = meanMatrix[2, 1 + i]
	sin, sout, diff, pvals[i] = statisticalTest(meanMatrix[2], d2i[i], 'smin', ptype='mean')
	#pvalue = np.minimum(1., pvalue * datasources.shape[0])

	print('%24s\t%.3f\t%.3f\t%.3f\t%.4f\t%.5f' % (d, salone, sin, sout, diff, pvals[i] ))
fig.savefig('../figures/clustering/per_data_source/' + species + '_smin.png', dpi=500)

print('\n\nbest score insted of mean\n')
print('Fmax')
print('%24s\t%s\t%s\t%s\t%s' % ('DATASOURCE', 'alone', 'in', 'out', 'diff'))
for i,d in enumerate(datasources):

	falone = meanMatrix[0, 1 + i]
	fin, fout, diff, pvals[i] = statisticalTest(meanMatrix[0], d2i[i], 'fmax', ptype='max')

	print('%24s\t%.3f\t%.3f\t%.3f\t%.4f\t%.5f' % (d, falone, fin, fout, diff, pvals[i] ))

print('\n\nSmin')
print('%24s\t%s\t%s\t%s\t%s' % ('DATASOURCE', 'alone', 'in', 'out', 'diff'))
for i,d in enumerate(datasources):

	sin, sout, diff, pvals[i] = statisticalTest(meanMatrix[2], d2i[i], 'smin', ptype='min')
	salone = meanMatrix[2, 1 + i]
	print('%24s\t%.3f\t%.3f\t%.3f\t%.4f\t%.5f' % (d, salone, sin, sout, diff, pvals[i] ))




plt.close('all')
