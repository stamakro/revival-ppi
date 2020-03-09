import pickle
from routinesgo import *
import numpy as np
import sys
from scipy.sparse import csr_matrix, save_npz

species = sys.argv[1]

try:
	ontology = sys.argv[2]
except IndexError:
	ontology ='P'


roots = {'P': 'GO:0008150', 'F': 'GO:0003674', 'C': 'GO:0005575'}

ontologyRoot = roots[ontology]

experimentalECs = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP']
highthroughputECs = ['HTP', 'HDA', 'HMP', 'HGI', 'HEP']
phylogeneticECs = ['IBA', 'IBD', 'IKR', 'IRD']
computationalECs = ['ISS', 'ISO', 'ISA', 'ISM', 'IGC', 'RCA']
authorECs = ['TAS', 'NAS']
curatorECs = ['IC']
electronicECs = ['IEA']

unwantedECs = electronicECs + computationalECs + ['NAS']
try:
	removeIpi = int(sys.argv[3])
except IndexError:
	removeIpi = 0

if removeIpi:
	unwantedECs += ['IPI']


if species != 'swissprot':
	path = '../data/' + species + '/annotations/'
	swissprot = set()
	with open('../data/' + species + '/sequences/sequences.list') as f:
		for line in f:
			swissprot.add(line[:-1])

else:
	path = '../blast/annotations/'
	swissprot = set()
	with open('../blast/swissprot.list') as f:
		for line in f:
			swissprot.add(line[:-1])



if removeIpi:
	prefix = ontology + '_' + 'noIPI_'
else:
	prefix = ontology + '_'

[allGenes, dag, mapping] = read_go_per_gene_ec(path + 'goa.gaf', obofile='../data/go/go-final.obo', excludedEcs=unwantedECs, allowedExplanations=dict())

[Y, termNames, geneNames] =  getLabelMatrix(allGenes, dag, ontology)

rootInd = termNames.index(ontologyRoot)
Y = np.delete(Y, rootInd, axis=1)

del termNames[rootInd]

if species != 'tomato':
	notReviewed = [i for i,g in enumerate(geneNames) if g not in swissprot]

	Y = np.delete(Y, notReviewed, axis=0)
	geneNames = np.delete(np.array(geneNames), notReviewed)


	emptyP = np.where(np.sum(Y, 1) == 0)[0]
	Y = np.delete(Y, emptyP, axis=0)

	geneNames = np.delete(geneNames, emptyP)

	emptyT = np.where(np.sum(Y, 0) == 0)[0]
	Y = np.delete(Y, emptyT, axis=1)

	termNames = np.delete(np.array(termNames), emptyT)

	#make sure there are no empty proteins
	assert np.min(np.sum(Y, 1)) > 0


else:
	#fix ids, to do
	gg2sol = dict()
	with open(path + 'gene2Solyc.map') as f:
		for line in f:

			fields = line.split('\n')[0].split('\t')

			if fields[0] in geneNames:
				if 'Solyc' in fields[1]:
					gg2sol[fields[0]] = fields[1]

				else:
					f2 = fields[2].split('|')
					for ff in f2:
						if 'Solyc' in ff:
							gg2sol[fields[0]] = ff
							break



assert Y.shape == (len(geneNames), len(termNames))


with open(path + prefix + 'geneNames.pkl', 'wb') as f:
	pickle.dump(geneNames, f)

with open(path + prefix + 'termNames.pkl', 'wb') as f:
	pickle.dump(termNames, f)

save_npz(path + prefix + 'Y.npz', csr_matrix(Y))
