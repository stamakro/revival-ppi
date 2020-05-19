import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pickle
import sys
from scipy.special import binom
from itertools import combinations
from sklearn.metrics import f1_score
from scipy.stats import spearmanr, pearsonr, wilcoxon
from routinesnetworks import *

def getFmaxPerProtein(method, species, nFolds=5, gbaIndex=-1, node2vecNet='biogrid'):
    path = '../experiments/' + species + '_P_median/'

    with open('../data/' + species + '/annotations/P_geneNames.pkl', 'rb') as f:
        geneNames = pickle.load(f)


    thresholds = np.linspace(0, 1., 51)
    N = geneNames.shape[0]
    gf = np.zeros((N,))

    predicted = np.zeros((N,), int)

    start = 0
    if method == 'blast':
        for fold in range(nFolds):
            print('fold %d' % fold)
            fullPath = path + 'fold' + str(fold) + '/'
            testNames = set()

            with open(fullPath + 'test.names') as f:
                for line in f:
                    testNames.add(line[:-1])

            with open(fullPath + 'blast.pkl', 'rb') as f:
                blastResults = pickle.load(f)

            termNames = blastResults['terms']
            del blastResults['terms']
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
            predicted[start + pred] = 1


            fmeasures = np.zeros((thresholds.shape[0], len(testNames)))

            for j, thr in enumerate(thresholds):

                if thr > np.max(Ypost):
                    break
                fmeasures[j] = f1_score(Ytrue.T, (Ypost.T >= thr).astype(int), average=None)


            best = np.argmax(np.mean(fmeasures, 1))
            gf[start:start+len(testNames)] = fmeasures[best]

            start += len(testNames)

    elif method == 'gba_1hop':
        for fold in range(nFolds):
            print('fold %d' % fold)
            fullPath = path + 'fold' + str(fold) + '/'

            testNames = set()
            with open(fullPath + 'test.names') as f:
                for line in f:
                    testNames.add(line[:-1])

            with open(fullPath + 'ppi/' + str(gbaIndex) + '_knn.pkl', 'rb') as f:
                dd = pickle.load(f)


            if gbaIndex > 0:
                Ypost = dd[0].toarray()
                with open(fullPath + 'ppi/0_knn.pkl', 'rb') as f:
                    Ytrue = pickle.load(f)['gt'].toarray()

            else:
                Ytrue = dd['gt'].toarray()
                Ypost = dd[str(gbaIndex)][0].toarray()

            Ypost[np.isnan(Ypost)] = 0

            fmeasures = np.zeros((thresholds.shape[0], len(testNames)))

            for j, thr in enumerate(thresholds):


                if thr > np.max(Ypost):
                    break
                fmeasures[j] = f1_score(Ytrue.T, (Ypost.T >= thr).astype(int), average=None)

            best = np.argmax(np.mean(fmeasures, 1))
            print(best, np.mean(fmeasures[best]))
            gf[start:start+len(testNames)] = fmeasures[best]

            pred = np.where(np.max(Ypost, axis=1) > 0.)[0]
            predicted[start + pred] = 1


            start += len(testNames)



    elif method == 'node2vec':
        for fold in range(nFolds):
            print('fold %d' % fold)
            fullPath = path + 'fold' + str(fold) + '/'

            testNames = set()
            with open(fullPath + 'test.names') as f:
                for line in f:
                    testNames.add(line[:-1])

            with open(fullPath + 'node2vec/' + node2vecNet + '.pkl', 'rb') as f:
                Ypost = pickle.load(f)


            with open(fullPath + 'ppi/0_knn.pkl', 'rb') as f:
                Ytrue = pickle.load(f)['gt'].toarray()

            Ypost[np.isnan(Ypost)] = 0

            fmeasures = np.zeros((thresholds.shape[0], len(testNames)))

            for j, thr in enumerate(thresholds):


                if thr > np.max(Ypost):
                    break
                fmeasures[j] = f1_score(Ytrue.T, (Ypost.T >= thr).astype(int), average=None)

            best = np.argmax(np.mean(fmeasures, 1))
            print(np.mean(fmeasures[best]))
            gf[start:start+len(testNames)] = fmeasures[best]

            pred = np.where(np.max(Ypost, axis=1) > 0.)[0]
            predicted[start + pred] = 1


            start += len(testNames)



    return gf, predicted


def getDegree(species):

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
    cv = KFold(n_splits=n_folds, shuffle=True, random_state=656391)

    np.random.seed(1901273)

    current = 0
    for fold, (train, test) in enumerate(cv.split(Aexp)):

        N = test.shape[0]

        degree[current:current+N] = np.sum(Aexp[test][:, train], 1)
        current += N



    return degree

if len (sys.argv) == 4:
    option1 = int(sys.argv[1])
    option2 = int(sys.argv[2])
    interactive = int(sys.argv[3])
    option3 = -1
elif len(sys.argv) == 5:
    option1 = int(sys.argv[1])
    option2 = int(sys.argv[2])
    option3 = int(sys.argv[3])
    interactive = int(sys.argv[4])


optionsList = dict()
optionsList[0] = ('yeast', 'gba_1hop', 0, None, 'PPI EXP')
optionsList[1] = ('yeast', 'gba_1hop', 5, None, 'PPI EXP + TM')
optionsList[2] = ('yeast', 'gba_1hop', 9, None, 'PPI EXP + Hom')
optionsList[3] = ('yeast', 'gba_1hop', 39, None, 'PPI EXP + TM + Hom') #best of string
optionsList[4] = ('yeast', 'blast', None, None, 'BLAST')
optionsList[5] = ('ecoli', 'gba_1hop', 129, None, 'PPI EXP + Fus + TMt + Hom')  #best of string
optionsList[6] = ('ecoli', 'gba_1hop', 0, None, 'PPI EXP')
optionsList[7] = ('arabidopsis', 'gba_1hop', 0, None, 'PPI EXP')
optionsList[8] = ('tomato', 'gba_1hop', 0, None, 'PPI EXP')
optionsList[9] = ('tomato', 'gba_1hop', 242, None, 'PPI EXP + Coext + EXPt + TM + Fus + TMt + Hom') #best of string
optionsList[10] = ('tomato', 'gba_1hop', 88, None, 'PPI EXP + TM + TMt + Hom') #equally good
optionsList[11] = ('arabidopsis', 'gba_1hop', 252, None, 'PPI EXP + TM + Cooc + Fus + Hom') #best of string
optionsList[12] = ('arabidopsis', 'gba_1hop', 39, None, 'PPI EXP + TM + Hom') #equally good
optionsList[13] = ('ecoli', 'blast', None, None, 'BLAST')
optionsList[14] = ('arabidopsis', 'blast', None, None, 'BLAST')
optionsList[15] = ('tomato', 'blast', None, None, 'BLAST')
optionsList[16] = ('yeast', 'node2vec', None, 'biogrid', 'PPI EXP n2v')
optionsList[17] = ('yeast', 'node2vec', None, 'biogrid_string', 'PPI EXP + TM + Hom n2v')
optionsList[18] = ('ecoli', 'node2vec', None, 'biogrid', 'PPI EXP n2v')
optionsList[19] = ('tomato', 'node2vec', None, 'biogrid', 'PPI EXP n2v')
optionsList[20] = ('arabidopsis', 'node2vec', None, 'biogrid', 'PPI EXP n2v')



species1, method1, gbaIndex1, node2vecNet1, name1 = optionsList[option1]
species2, method2, gbaIndex2, node2vecNet2, name2 = optionsList[option2]
assert species1 == species2
species = species1
del species1, species2

if option3 != -1:
    species3, method3, gbaIndex3, node2vecNet3, name3 = optionsList[option3]
    assert species3 == species
    del species3

#textmining:5

nFolds = 5

if not interactive:
    f1, pred1 = getFmaxPerProtein(method1, species, gbaIndex=gbaIndex1, node2vecNet=node2vecNet1)
    f2, pred2 = getFmaxPerProtein(method2, species, gbaIndex=gbaIndex2, node2vecNet=node2vecNet2)
    if option3 != -1:
        f3, pred3 = getFmaxPerProtein(method3, species, gbaIndex=gbaIndex3, node2vecNet=node2vecNet3)



m = np.minimum(np.min(f1), np.min(f2))
M = np.maximum(np.max(f1), np.max(f2))

degrees = getDegree(species)

r_p1, _ = pearsonr(degrees, f1)
r_s1, _ = spearmanr(degrees, f1)

r_p2, _ = pearsonr(degrees, f2)
r_s2, _ = spearmanr(degrees, f2)

space = max(len(name1), len(name2))

print('Correlation with degree:')
print('%0*s\tPearson\t\tSpearman' % (space, 'Model'))
print('%0*s\t%.3f\t\t%.3f' % (space, name1, r_p1, r_s1))
print('%0*s\t%.3f\t\t%.3f' % (space, name2, r_p2, r_s2))



from scipy.stats import gaussian_kde, linregress
a, b, r, p, stderr = linregress(f1, f2)

xmin = m-0.05
xmax = 1.05*M
ymin = m-0.05
ymax = 1.05*M

fig = plt.figure()
ax = fig.add_subplot(111)

dd = np.vstack((f1, f2))
X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
positions = np.vstack([X.ravel(), Y.ravel()])
kernel = gaussian_kde(dd)
Z = np.reshape(kernel(positions).T, X.shape)
fob = ax.pcolormesh(X, Y, Z, cmap='Oranges')
xx = np.linspace(0.0, 1.0)
ax.plot(xx, xx, 'k--')
plt.colorbar(fob)



ax.scatter(f1, f2, color='C0', alpha=0.1, s=12, edgecolor='k')


ax.plot(xx, a*xx + b, 'r')

ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_xlabel(name1)
ax.set_ylabel(name2)
ax.set_title('Fmax ' + species)
ax.set_aspect('equal')


plt.tight_layout()
fig.savefig('../figures/others/' + species + '_' + name1 + '_' + name2 + '.png', dpi=300)

'''
fig = plt.figure()
ax = fig.add_subplot(1,2,1)
ax.scatter(np.log(degrees+1), f1, edgecolor='k')
ax.set_xlabel('Node degree')
ax.set_ylabel('Fmax ' + name1)
#plt.xscale('log')


ax = fig.add_subplot(1,2,2)
ax.scatter(np.log(degrees+1), f2, edgecolor='k')
ax.set_xlabel('Node degree')
ax.set_ylabel('Fmax ' + name2)
#plt.xscale('log')

plt.tight_layout()
fig.savefig('../figures/others/' + species + '_degree_gba_n2v.png')
'''

fig = plt.figure()
if species == 'tomato':
    bins = [-0.5, 0.5, 12.5]
    labels = ['0', '1-12']

elif species == 'yeast':
    bins = [-0.5, 0.5, 2.5, 5.5, 10.5, 20.5, 50.5, 100.5, 200.5]
    labels = ['0', '1-2', '3-5', '6-10', '11-20', '21-50', '51-100', '101-200', '>200']

elif species == 'arabidopsis':
    bins = [-0.5, 0.5, 2.5, 5.5, 10.5]
    labels = ['0', '1-2', '3-5', '6-10', '>10' ]

elif species == 'ecoli':
    bins = [-0.5, 0.5, 2.5, 5.5, 10.5]
    labels = ['0', '1-2', '3-5', '6-10', '>10' ]

indices = np.digitize(degrees, bins)
indexVals = np.unique(indices)

if option3 == -1:
    x1 = 3 * np.arange(len(bins)-1) - 0.5
    x2 = 3 * np.arange(len(bins)-1) + 0.5

    ff1 = []
    ff2 = []

    for ii in indexVals:
        d1 = f1[np.where(indices == ii)]
        d2 = f2[np.where(indices == ii)]
        ff1.append(d1)
        ff2.append(d2)
        print(wilcoxon(d1, d2))

    ax = fig.add_subplot(111)
    bp1 = ax.boxplot(ff1, whis=[5, 95], positions=x1, patch_artist=True, medianprops={'color':'k'}, showfliers=False)
    bp2 = ax.boxplot(ff2, whis=[5, 95], positions=x2, patch_artist=True, medianprops={'color':'k'}, showfliers=False)

    for patch in bp1['boxes']:
        patch.set_facecolor('C0')
        #patch.set_alpha(0.5)

    for patch in bp2['boxes']:
        patch.set_facecolor('C1')
        #patch.set_alpha(0.5)

    ax.set_xticks(3 * np.arange(len(bins) - 1))
    ax.set_xticklabels(['0', '1-2', '3-5', '6-10', '11-20', '>20'])

    ax.set_xlim(-1, np.max(x2)+0.5)

    ax.scatter([], [], color='C0', marker='s', s=30, edgecolor='k', label=name1)
    ax.scatter([], [], color='C1', marker='s', s=30, edgecolor='k', label=name2)
    ax.set_ylim(-0.2, 1.05)

    ax.legend(loc='lower right')
    fig.savefig('../figures/others/' + species + '_degree_' + name1 + '_' + name2 + '_perlevel.png')

else:
        NN = 4.0
        x1 = NN * np.arange(len(labels)) - 0.5
        x2 = NN * np.arange(len(labels)) + 0.5
        x3 = NN * np.arange(len(labels)) + 1.5


        ff1 = []
        ff2 = []
        ff3 = []

        for ii in indexVals:
            d1 = f1[np.where(indices == ii)]
            d2 = f2[np.where(indices == ii)]
            d3 = f3[np.where(indices == ii)]

            ff1.append(d1)
            ff2.append(d2)
            ff3.append(d3)
            print(wilcoxon(d1, d2))

        ax = fig.add_subplot(111)
        bp1 = ax.boxplot(ff1, whis=[5, 95], positions=x1, patch_artist=True, medianprops={'color':'k'}, showfliers=False)
        bp2 = ax.boxplot(ff2, whis=[5, 95], positions=x2, patch_artist=True, medianprops={'color':'k'}, showfliers=False)
        bp3 = ax.boxplot(ff3, whis=[5, 95], positions=x3, patch_artist=True, medianprops={'color':'k'}, showfliers=False)

        for patch in bp1['boxes']:
            patch.set_facecolor('C0')
            #patch.set_alpha(0.5)

        for patch in bp2['boxes']:
            patch.set_facecolor('C1')
            #patch.set_alpha(0.5)

        for patch in bp3['boxes']:
            patch.set_facecolor('c')
            #patch.set_alpha(0.5)

        #ax.set_xticks(3.5 * np.arange(len(bins) - 1))
        ax.set_xticks(x2)
        ax.set_xticklabels(labels)

        ax.set_xlim(-1, np.max(x3)+0.5)

        for jj, ss in zip(x2, ff1):
            ax.text(jj-0.25, 1.05, str(len(ss)))
  
        ax.scatter([], [], color='C0', marker='s', s=30, edgecolor='k', label=name1)
        ax.scatter([], [], color='C1', marker='s', s=30, edgecolor='k', label=name2)
        ax.scatter([], [], color='c', marker='s', s=30, edgecolor='k', label=name3)
        ax.set_ylim(-0.3, 1.1)

        ax.legend(loc='lower right')
        ax.set_xlabel('Number of annotated neighbors in EXP network')
        ax.set_ylabel('F1 score')

        fig.savefig('../figures/others/' + species + '_degree_' + name1 + '_' + name2 + '_' + name3 + '_perlevel.png')
        fig.savefig('../figures/others/' + species + '_degree_' + name1 + '_' + name2 + '_' + name3 + '_perlevel.eps', dpi=1200)
