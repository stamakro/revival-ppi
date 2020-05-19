import pickle
import sys
from routinesnetworks import getPPInetwork
import numpy as np

np.random.seed(17042010)
species = sys.argv[1]

try:
    annoPrefix = sys.argv[2]
except IndexError:
    annoPrefix = 'P_'


with open('../data/' + species + '/annotations/' + annoPrefix + 'geneNames.pkl', 'rb') as f:
    geneNames = pickle.load(f)


seqDict = dict()
save = False

with open('../data/' + species + '/sequences/sequences.fasta') as f:
    for line in f:
        if line[0] == '>':
            if save:
                seqDict[protein] = seq
            protein = line.split('|')[1]
            assert protein not in seqDict

            save = False
            seq = ''

            if protein in geneNames:
                save = True

        else:
            seq += line[:-1]

if save:
    seqDict[protein] = seq

assert len(seqDict) == len(geneNames)


Aexp = getPPInetwork(species, 'biogrid')
assert np.max(np.abs((Aexp - Aexp.T)))  < 1e-10

assert Aexp.shape[0] == len(seqDict)

ind = np.triu_indices_from(Aexp)

pos = np.where(Aexp[ind])[0]
neg = np.where(Aexp[ind] == 0)[0]
neg = np.random.permutation(neg)[:pos.shape[0]]

p1 = np.unique(ind[0][pos])
p2 = np.unique(ind[1][pos])
n1 = np.unique(ind[0][neg])
n2 = np.unique(ind[1][neg])

#assert np.union1d(np.union1d(np.union1d(p1,p2), n1), n2).shape[0] == len(geneNames)

with open('seq_ppi/' + species + '/biogrid/protein.dictionary.training.tsv', 'w') as f:
    for g in geneNames:
        f.write(g + '\t' + seqDict[g] + '\n')

with open('seq_ppi/' + species + '/biogrid/protein.actions.training.tsv', 'w') as f:
    for i in range(pos.shape[0]):
        p1 = ind[0][pos[i]]
        p2 = ind[1][pos[i]]
        f.write(geneNames[p1] + '\t' + geneNames[p2] + '\t1\n')

        p1 = ind[0][neg[i]]
        p2 = ind[1][neg[i]]
        f.write(geneNames[p1] + '\t' + geneNames[p2] + '\t0\n')
