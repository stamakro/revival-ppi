import pickle
import sys
import os

expPath = sys.argv[1]
nFolds = 5

testSets = [set() for _ in range(nFolds)]

for fold in range(nFolds):

    with open(expPath + '/fold' + str(fold) + '/test.names') as f:
        for line in f:
            testSets[fold].add(line[:-1])

annoType = expPath.split('_')
species = annoType[0].split('/')[-1]

#if len(annoType) != 3:
#    raise NotImplementedError
if 'noIPI' in annoType:
    annoPrefix = annoType[1] + '_' + annoType[2] + '_'
else:
    annoPrefix = annoType[1] + '_'

with open('../blast/annotations/' + annoPrefix + 'geneNames.pkl', 'rb') as f:
    annoGenes = set(pickle.load(f))

fileDsTrn = [open('../blast/subsets/' + species + 'training_fold' + str(i) + '.fasta', 'w') for i in range(nFolds)]
fileDsTst = [open('../blast/subsets/' + species + 'test_fold' + str(i) + '.fasta', 'w') for i in range(nFolds)]

with open('../blast/swissprot.fasta') as f:
    for line in f:
        if line[0] == '>':
            proteinId = line.split('|')[1]


            if proteinId not in annoGenes:
                writeTrn = [False for _ in range(nFolds)]
                writeTst = [False for _ in range(nFolds)]
            else:
                writeTrn = [True for _ in range(nFolds)]
                writeTst = [False for _ in range(nFolds)]

                for i in range(nFolds):
                    if proteinId in testSets[i]:
                        writeTrn[i] = False
                        writeTst[i] = True


        for i, (fwtr, fwts) in enumerate(zip(fileDsTrn, fileDsTst)):
            if writeTrn[i]:
                fwtr.write(line)

            if writeTst[i]:
                fwts.write(line)


for ff in fileDsTrn:
    ff.close()
for ff in fileDsTst:
    ff.close()
