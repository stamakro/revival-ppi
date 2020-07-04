import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
import pickle


#figure 1, bar chart

species = ['S. cerevisiae', 'E. coli', 'A. thaliana', 'S. lycopersicum']
names = ['EXP', 'EXP +\n BLAST', 'EXP +\n STRING', 'EXP +\n STRING +\n BLAST', 'EXP +\n SEQ']
y = np.array([
[0.42, 0.50, 0.49, 0.54, 0.33],
[0.25, 0.45, 0.46, 0.58, 0.32],
[0.19, 0.43, 0.48, 0.54, 0.29],
[0.08, 0.35, 0.61, 0.60, 0.35]
])

yerr = np.array([
[0.007, 0.007, 0.005, 0.006, 0.005],
[0.011, 0.020, 0.008, 0.021, 0.007],
[0.007, 0.006, 0.004, 0.008, 0.004],
[0.007, 0.025, 0.045, 0.047, 0.016]
])


naive = np.array([0.31, 0.29, 0.28, 0.23])
naiveErr = np.array([0.004, 0.009, 0.005, 0.019])

blast = np.array([0.35, 0.43, 0.42, 0.34])
blastErr = np.array([0.004, 0.021, 0.007, 0.025])


x = np.arange(1, 1+y.shape[1])

y2 = np.array([
[0.50, 0.59],
[0.28, 0.50],
[0.23, 0.50],
[0.08, 0.61]
])

y2err = np.array([
[0.012, 0.009],
[0.011, 0.012],
[0.008, 0.005],
[0.042, 0.042]
])

x = 1.5 * x
#i = 3
letters = ['a', 'c', 'b', 'd']
locations = [1, 2, 4, 5]
fig = plt.figure(figsize=(12,9))
lw = 0.8

for i in range(len(species)):

    ax = fig.add_subplot(2,3,locations[i])

    ax.bar(x, y[i], yerr=yerr[i], color=['C0', 'C7', 'C9', 'C8', 'g'], edgecolor='k')

    #ax.bar(x[-3:-1], y2[i] - y[i,-3:-1], bottom=y[i, -3:-1], yerr=y2err[i], color=['C0', 'C1'], edgecolor='k', hatch='\\',  label=['GBA' 'node2vec'])
    for j, hh in enumerate(y2[i]):
        if hh >= y[i, 2*j]:
            ax.bar(x[2*j], hh - y[i,2*j], bottom=y[i,2*j], yerr=y2err[i,j], color='C1', hatch='\\', edgecolor='k')
        else:
            ax.bar(x[2*j], hh, bottom=0, yerr=y2err[i,j], color='C1', hatch='\\', edgecolor='k')


    ax.axhline(naive[i], color='k')
    ax.axhline(naive[i]+naiveErr[i], color='k', linestyle='--', linewidth=lw)
    ax.axhline(naive[i]-naiveErr[i], color='k', linestyle='--', linewidth=lw)

    ax.axhline(blast[i], color='C3')
    ax.axhline(blast[i]+blastErr[i], color='C3', linestyle='--', linewidth=lw)
    ax.axhline(blast[i]-blastErr[i], color='C3', linestyle='--', linewidth=lw)

    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_title(letters[i] + ') ' + species[i])
    ax.set_ylabel('protein-centric Fmax')



ax = fig.add_subplot(2, 3, 3)
ee = []
for species in ['yeast', 'ecoli', 'arabidopsis', 'tomato']:
    with open('../data/' + species + '/degrees/0.pkl', 'rb') as f:
        d = pickle.load(f)
        E = np.sum(d/2)
        N = d.shape[0] * (d.shape[0] - 1) / 2
        ee.append(100*E/N)

perc = (y[:,2] - y[:, 0]) / y[:, 0]
pa = 1 - ((1 - y[:,2]) / (1 - y[:,0]))

ax.scatter(ee, pa, s=60, edgecolor='k')

slope, intercept, rho, pvalue, stderr = linregress(ee, pa)
xx = np.linspace(np.min(ee), np.max(ee), 3)
ax.plot(xx, intercept + slope * xx, color='r')

ax.set_xlabel('Percentage of proteins known to interact')
ax.set_ylabel('Prediction Advantage in Fmax')

ax.set_title('e) PA of GBA classifier')


ax = fig.add_subplot(2, 3, 6)

Nproteins = [4997, 2869, 10648, 651]
speciesNames = ['S. cerevisiae', 'E. coli', 'A. thaliana', 'S. lycopersicum']
methods = ['BLAST', 'EXP GBA', 'EXP node2vec', 'EXP + STRING GBA', 'EXP + SEQ GBA']


coverage = np.array([
[0.877, 0.994, 0.996, 0.998, 1.],
[0.852, 0.749, 0.769, 0.974, 0.999],
[0.973, 0.532, 0.566, 0.958, 0.975],
[0.989, 0.029, 0.031, 0.796, 0.977]
])

coverageStd = np.array([
[0.007, 0.002, 0.002, 0.001, 0.],
[0.013, 0.013, 0.016, 0.005, 0.001],
[0.003, 0.017, 0.011, 0.002, 0.002],
[0.007, 0.015, 0.017, 0.026, 0.019]
])

w = 0.8

x = np.arange(len(speciesNames)) * len(methods)
offset = np.arange(len(methods)) * w
#colors = ['C' + str(i) for i in range(len(methods))]
colors = ['C3','C0', 'C1', 'c', 'g']


for j in range(coverage.shape[1]):
    ax.bar(x + offset[j], coverage[:,j], yerr=coverageStd[:,j], color=colors[j], label=methods[j], edgecolor='k')


centers = np.median(offset) + x
ax.set_xticks(centers)
ax.set_xticklabels(speciesNames, rotation=45)

for (xx, s) in zip(x, Nproteins):
    ax.text(xx, 1.05, str(s))

ax.set_ylim(0., 1.1)
#ax.legend(loc='upper right')
#ax.legend(loc='best', bbox_to_anchor=(0.0, 1.5, 1.0, 0.4))
ax.set_ylabel('Fraction of protein annotated')

ax.set_title('f) Coverage')

plt.tight_layout()
fig.savefig('../figures/main-paper/bar.png')
fig.savefig('../figures/main-paper/bar.eps', dpi=1200)


species = ['yeast', 'ecoli', 'arabidopsis', 'tomato']
Nbins = [30, 30, 50, 15]
speciesNames = ['(a) S. cerevisiae', '(b) E. coli', '(c) A. thaliana', '(d) S. lycopersicum']

fig = plt.figure()

bb = np.array([-0.5, 0.5, 2.5, 5.5, 10.5, 20.5, 50.5, 100.5, 200.5, 500.5, 1000.5, 2000.5])

for i, (s, sn) in enumerate(zip(species, speciesNames),1):
	ax = fig.add_subplot(2,2,i)

	with open('../data/' + s + '/degrees/0.pkl', 'rb') as f:
		dd = pickle.load(f)


	print(sn)
	print(np.percentile(dd, [1, 5, 10, 25, 50, 75, 90, 95, 99]))

	#ax.hist(np.log10(dd+1), bins=30, log=True, color='C0', edgecolor='k', density=True)
	ax.hist(dd, bins=Nbins[i-1], log=True, color='C0', edgecolor='k', density=True)


	ax.set_title(sn)
	ax.set_xlabel('Degree')
	ax.set_ylabel('Log probability density')


	#ticks = []


plt.tight_layout()

fig.savefig('../figures/main-paper/degrees.png')

fig.savefig('../figures/main-paper/degrees.eps', dpi=1200)





plt.close('all')
