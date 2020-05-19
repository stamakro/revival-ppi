import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pickle
import sys
from scipy.special import binom

experimentPath = sys.argv[1]
stuff = experimentPath.split('/')[2]
stuff2 = stuff.split('_')
species = stuff2[0]

fignames = ['coverage', 'global', 'local']

classifiers = ['gba_1hop', 'gba+blast_1hop', 'n2v_5nn']
classifiers = ['gba_1hop']
n_folds = 5

if species == 'tomato':
	ds = 8
else:
	ds = 9

metrics = ['coverage', 'global_fmax', 'local_fmax', 'global_smin', 'local_smin']

mS = np.zeros((ds+1,))
mF = np.zeros((ds+1,))
mS_local = np.zeros((ds+1,))
mF_local = np.zeros((ds+1,))
cov = np.zeros((ds+1,))


with open(experimentPath + 'naive.pkl', 'rb') as f:
	naive = pickle.load(f)


for clf in classifiers:
    for fold in range(n_folds):
        with open(experimentPath + 'fold' + str(fold) + '/performance_' + clf + '.pkl', 'rb') as f:
            data = pickle.load(f)

        if fold == 0:
            bigMatrix = np.zeros((n_folds, len(metrics), data[0].shape[0]))

        for i, v in enumerate(data):
            bigMatrix[fold, i] = v

    meanMatrix = np.mean(bigMatrix, axis=0)
    stdMatrix = np.std(bigMatrix, axis=0, ddof=1)

    if clf == 'gba_1hop' or clf == 'gba_2hop' :
        blast = meanMatrix[:,0]
        biogrid = meanMatrix[:,1]
        blastAndBiogrid = meanMatrix[:,2]
        stringStart = 3

    elif clf == 'gba+blast_1hop':
        blast = meanMatrix[:,0]
        biogrid = meanMatrix[:,2]
        blastAndBiogrid = meanMatrix[:,2]
        stringStart = 3


    else:
        biogrid = meanMatrix[:,0]
        stringStart = 1

    figures = [plt.figure() for _ in range(3)]
    #0 coverage, 1 global, 2 local

    axCov = figures[0].add_subplot(1, 1, 1)
    ax_g_f = figures[1].add_subplot(1, 2, 1)
    ax_g_s = figures[1].add_subplot(1, 2, 2)

    ax_l_f = figures[2].add_subplot(1, 2, 1)
    ax_l_s = figures[2].add_subplot(1, 2, 2)


    counter = stringStart

    for i in range(ds+1):
        if i == 0:
            ax_g_f.scatter(0, biogrid[1], color='C0', edgecolor='k')
            ax_g_f.axvline(1, color='k', linestyle=':', alpha=0.3)

            ax_g_s.scatter(0, biogrid[3], color='C0', edgecolor='k')
            ax_g_s.axvline(1, color='k', linestyle=':', alpha=0.3)
            mF[0] = biogrid[1]
            mS[0] = biogrid[3]

            ax_l_f.scatter(0, biogrid[2], color='C0', edgecolor='k')
            ax_l_f.axvline(1, color='k', linestyle=':', alpha=0.3)

            ax_l_s.scatter(0, biogrid[4], color='C0', edgecolor='k')
            ax_l_s.axvline(1, color='k', linestyle=':', alpha=0.3)
            mF_local[0] = biogrid[2]
            mS_local[0] = biogrid[4]

            axCov.scatter(0, biogrid[0], color='C0', edgecolor='k')
            axCov.axvline(1, color='k', linestyle=':', alpha=0.3)
            cov[0] = biogrid[0]

        else:
            iterations = int(binom(ds, i))

            mF[i] = np.mean(meanMatrix[1, counter:counter+iterations])
            mS[i] = np.mean(meanMatrix[3, counter:counter+iterations])
            mF_local[i] = np.mean(meanMatrix[2, counter:counter+iterations])
            mS_local[i] = np.mean(meanMatrix[4, counter:counter+iterations])
            cov[i] = np.mean(meanMatrix[0, counter:counter+iterations])

            ax_g_f.scatter(2*i + 1.5 * np.random.rand(iterations) - 0.75, meanMatrix[1, counter:counter+iterations], color='C0', edgecolor='k' )
            ax_g_f.axvline(2*i+1, color='k', linestyle=':', alpha=0.3)

            ax_g_s.scatter(2*i + 1.5 * np.random.rand(iterations) - 0.75, meanMatrix[3, counter:counter+iterations], color='C0', edgecolor='k' )
            ax_g_s.axvline(2*i+1, color='k', linestyle=':', alpha=0.3)

            ax_l_f.scatter(2*i + 1.5 * np.random.rand(iterations) - 0.75, meanMatrix[2, counter:counter+iterations], color='C0', edgecolor='k' )
            ax_l_f.axvline(2*i+1, color='k', linestyle=':', alpha=0.3)

            ax_l_s.scatter(2*i + 1.5 * np.random.rand(iterations) - 0.75, meanMatrix[4, counter:counter+iterations], color='C0', edgecolor='k' )
            ax_l_s.axvline(2*i+1, color='k', linestyle=':', alpha=0.3)

            axCov.scatter(2*i + 1.5 * np.random.rand(iterations) - 0.75, meanMatrix[0, counter:counter+iterations], color='C0', edgecolor='k' )
            axCov.axvline(2*i+1, color='k', linestyle=':', alpha=0.3)

            counter += iterations


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
    #ax2.set_ylim(0, np.max(meanMatrix[3])*1.15)
    ax_g_f.axhline(naive[0][0], color='k', label='naive')
    ax_g_s.axhline(naive[1][0], color='k', label='naive')
    ax_g_f.axhline(blast[1], color='C3', label='blast')
    ax_g_s.axhline(blast[3], color='C3', label='blast')
    ax_g_f.axhline(blastAndBiogrid[1], color='C2', linestyle='--', label='blast+biogrid')
    ax_g_s.axhline(blastAndBiogrid[3], color='C2', linestyle='--', label='blast+biogrid')

    ax_g_f.legend()
    ax_g_s.legend()
    #figures[1].suptitle(species + ', ' + clf)
    figures[1].tight_layout()


    ax_l_f.plot(np.arange(0,2*mF_local.shape[0], 2), mF_local, color='C1')
    ax_l_s.plot(np.arange(0, 2*mS_local.shape[0], 2), mS_local, color='C1')
    ax_l_f.set_xticks([2*i for i in range(ds+2)])
    ax_l_f.set_xticklabels([str(i) for i in range(ds+2)])
    ax_l_s.set_xticks([2*i for i in range(ds+2)])
    ax_l_s.set_xticklabels([str(i) for i in range(ds+2)])
    ax_l_f.set_xlabel('# of data sources', fontsize=20)
    ax_l_f.set_ylabel('local Fmax', fontsize=20)
    ax_l_s.set_xlabel('# of data sources', fontsize=20)
    ax_l_s.set_ylabel('local Smin', fontsize=20)
    ax_l_f.set_xlim(-1, 2*ds+2)
    ax_l_s.set_xlim(-1, 2*ds+2)

    if clf == 'gba_1hop' and species == 'yeast':
        ax_l_f.set_ylim(0.3, 0.6)
    #ax2.set_ylim(0, np.max(meanMatrix[3])*1.15)
    ax_l_f.axhline(naive[0][0], color='k', label='naive')
    ax_l_s.axhline(naive[1][0], color='k', label='naive')
    ax_l_f.axhline(blast[2], color='C3', label='blast')
    ax_l_s.axhline(blast[4], color='C3', label='blast')
    ax_l_f.axhline(blastAndBiogrid[2], color='C2', linestyle='--', label='blast+biogrid')
    ax_l_s.axhline(blastAndBiogrid[4], color='C2', linestyle='--', label='blast+biogrid')

    ax_l_f.legend()
    ax_l_s.legend()
    figures[2].suptitle(species + ', ' + clf)


    axCov.plot(np.arange(0,2*cov.shape[0], 2), cov, color='C1')
    axCov.set_xticks([2*i for i in range(ds+2)])
    axCov.set_xticklabels([str(i) for i in range(ds+2)])
    axCov.set_xlabel('# of data sources', fontsize=20)
    axCov.set_ylabel('coverage', fontsize=20)
    axCov.set_xlim(-1, 2*ds+2)
    axCov.set_xlim(-1, 2*ds+2)

    axCov.axhline(blast[0], color='C3', label='blast')
    axCov.axhline(blastAndBiogrid[0], color='C2', linestyle='--', label='blast+biogrid')

    axCov.legend()
    figures[0].suptitle(species + ', ' + clf)


    for jj, fig in enumerate(figures):
        fig.savefig('../figures/' + species + '_' + clf + '_' + fignames[jj] + '.png')
        fig.savefig('../figures/' + species + '_' + clf + '_' + fignames[jj] + '.eps', dpi=1200)


plt.close('all')
