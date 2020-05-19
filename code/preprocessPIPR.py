import pickle
import sys


species = sys.argv[1]

try:
    annoPrefix = sys.argv[2]
except IndexError:
    annoPrefix = 'P_'


with open('../data/' + species + '/annotations/' + annoPrefix + 'geneNames.pkl', 'rb') as f:
    geneNames = pickle.load(f)


seqDict = dict()
save = False

if species != 'tomato':
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
else:
    with open('../data/' + species + '/sequences/sequences_tomatoIDs.fasta') as f:
        for line in f:
            if line[0] == '>':
                if save:
                    seqDict[protein] = seq
                protein = line[1:].split(' ')[0].split('.')[0]
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


with open('seq_ppi/' + species + '/preprocessed/protein.dictionary.test.tsv', 'w') as f:
    for g in geneNames:
        f.write(g + '\t' + seqDict[g] + '\n')

with open('seq_ppi/' + species + '/preprocessed/protein.actions.test.tsv', 'w') as f:
    for i in range(len(geneNames) - 1):
        for j in range(i+1, len(geneNames)):
            f.write(geneNames[i] + '\t' + geneNames[j] + '\t1\n')
