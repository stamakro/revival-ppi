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

path = '../data/' + species + '/annotations/'

if removeIpi:
	prefix = ontology + '_' + 'noIPI_'
else:
	prefix = ontology + '_'

[allGenes, dag, mapping] = read_go_per_gene_ec(path + 'goa.gaf', obofile='../data/go/go-final.obo', excludedEcs=unwantedECs, allowedExplanations=dict())

[Y, termNames, geneNames] =  getLabelMatrix(allGenes, dag, ontology)

rootInd = termNames.index(ontologyRoot)
Y = np.delete(Y, rootInd, axis=1)

del termNames[rootInd]

gdel = np.where(np.sum(Y, 1) == 0)[0]
Y = np.delete(Y, gdel, axis=0)

for tt in sorted(gdel)[::-1]:
	del geneNames[tt]

with open(path + prefix + 'geneNames.pkl', 'wb') as f:
	pickle.dump(geneNames, f)

with open(path + prefix + 'termNames.pkl', 'wb') as f:
	pickle.dump(termNames, f)

save_npz(path + prefix + 'Y.npz', csr_matrix(Y))
