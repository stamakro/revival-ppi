import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



with open('../results/downsampling_yeast_random.pkl', 'rb') as f:
	stuffR = pickle.load(f)


with open('../results/downsampling_yeast_degree.pkl', 'rb') as f:
	stuffD = pickle.load(f)



fracMiss = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7 , 0.8 , 0.9, 0.99])


fmaxR = np.mean(stuffR['gf'], axis=2)
sminR = np.mean(stuffR['gs'], axis=2)
cR = np.mean(stuffR['c'], axis=2)

fmaxD = np.mean(stuffD['gf'], axis=2)
sminD = np.mean(stuffD['gs'], axis=2)
cD = np.mean(stuffD['c'], axis=2)


fmaxRmean = np.mean(fmaxR, axis=0)
fmaxRstd = np.std(fmaxR, axis=0, ddof=1)

sminRmean = np.mean(sminR, axis=0)
sminRstd = np.std(sminR, axis=0, ddof=1)

cRmean = np.mean(cR, axis=0)
cRstd = np.std(cR, axis=0, ddof=1)

fmaxDmean = np.mean(fmaxD, axis=0)
fmaxDstd = np.std(fmaxD, axis=0, ddof=1)

sminDmean = np.mean(sminD, axis=0)
sminDstd = np.std(sminD, axis=0, ddof=1)

cDmean = np.mean(cD, axis=0)
cDstd = np.std(cD, axis=0, ddof=1)



fig = plt.figure(figsize=(12, 7))
ax = fig.add_subplot(1,3,1)

ax.scatter(fracMiss, fmaxRmean, color='y',s=5)
ax.errorbar(fracMiss, fmaxRmean,yerr=fmaxRstd, color='y', ecolor='k', capsize=4, label='uniform')

ax.scatter(fracMiss, fmaxDmean, color='m', s=5)
ax.errorbar(fracMiss, fmaxDmean,yerr=fmaxDstd, color='m', ecolor='k', capsize=4, label='degree')

ax.set_xlabel('Fraction of missing edges\n(a)', fontsize=14)
ax.set_ylabel('Fmax', fontsize=14)
ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
plt.legend()

ax = fig.add_subplot(1,3,2)

ax.scatter(fracMiss, sminRmean, color='y',s=5)
ax.errorbar(fracMiss, sminRmean,yerr=sminRstd, color='y', ecolor='k', capsize=4, label='uniform')

ax.scatter(fracMiss, sminDmean, color='m', s=5)
ax.errorbar(fracMiss, sminDmean, yerr=sminDstd, color='m', ecolor='k', capsize=4, label='degree')

ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_xlabel('Fraction of missing edges\n(b)', fontsize=14)
ax.set_ylabel('Smin', fontsize=14)
plt.legend()

ax = fig.add_subplot(1,3,3)

ax.scatter(fracMiss, cRmean, color='y',s=5)
ax.errorbar(fracMiss, cRmean,yerr=cRstd, color='y', ecolor='k', capsize=4, label='uniform')

ax.scatter(fracMiss, cDmean, color='m', s=5)
ax.errorbar(fracMiss, cDmean, yerr=cDstd, color='m', ecolor='k', capsize=4, label='degree')

ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_xlabel('Fraction of missing edges\n(c)', fontsize=14)
ax.set_ylabel('Coverage', fontsize=14)
plt.legend()

fig.savefig('../figures/tmp_downsample.png',dpi=500)
fig.savefig('../figures/tmp_downsample.eps',dpi=1200)





