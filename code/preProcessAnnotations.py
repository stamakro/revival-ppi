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
	with open('../data/tomato/annotations/gene2solyc.pkl', 'rb') as f:
		g2sol = pickle.load(f)

	notReviewed = [i for i,g in enumerate(geneNames) if g not in g2sol]

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

	geneNames = np.array([g2sol[g] for g in geneNames])

	genes, counts = np.unique(geneNames, return_counts = True)
	tobedel = []
	for g, cc in zip(genes, counts):
		if cc > 1:
			assert cc == 2
			locations = np.where(geneNames == g)[0]
			assert (Y[locations[0]] == Y[locations[1]]).all()
			tobedel.append(locations[1])

	Y = np.delete(Y, tobedel, axis=0)
	geneNames = np.delete(geneNames, tobedel)
	gene2row = dict()

	for i, g in enumerate(geneNames):
	    gene2row[g] = i
	genes = set(geneNames)


	reviewed = dict()
	with open('../data/tomato/sequences/sequences3.tab') as f:
	    for i, line in enumerate(f):
	        if i == 0:
	            continue

	        fields = line.split('\t')
	        uni = fields[0]
	        plant = fields[-1]
	        status = fields[2]

	        assert uni not in reviewed

	        if status == 'reviewed':
	            reviewed[uni] = True
	        else:
	            reviewed[uni] = False


	g2u = dict()
	u2g = dict()
	foundHits = set()

	foundExact = set()

	with open('../data/tomato/sequences/temp-blast/blastout.txt') as f:
	    for line in f:
	        if line[0] != '#':
	            fields = line.split()

	            plant = fields[0].split('.')[0]
	            uni = fields[1].split('|')[1]
	            identity = float(fields[2])
	            evalue = float(fields[-2])
	            length = int(fields[3])
	            qstart = int(fields[6])
	            qend = int(fields[7])
	            qCov = length / (qend - qstart + 1)

	            #if plant in foundExact:
	            #    #assert identity < 100
	            #    continue


	            if identity > 99 and qCov > 0.95:
	                #if np.abs(qCov - 1) > 1e-6:
	                #    print(line, end='')

	                if identity == 100.:
	                    foundExact.add(plant)

	                foundHits.add(plant)
	                if plant not in g2u:
	                    g2u[plant] = {uni}
	                else:
	                    g2u[plant].add(uni)

	                if uni not in u2g:
	                    u2g[uni] = {plant}
	                else:
	                    u2g[uni].add(plant)

	g2uFinal = dict()
	u2gFinal = dict()

	tobedelG = set()
	tobedelU = set()

	for g in g2u:
	    if g in tobedelG or g in g2uFinal:
	        continue


	    if len(g2u[g]) == 1:
	        u = g2u[g].pop()
	        g2u[g] = {u}
	        if len(u2g[u]) == 1:
	            #for this case, solyc id matches uniquely to a uniprot id and vice versa.
	            #these are the easy, clean cases
	            gg = u2g[u].pop()
	            u2g[u] = {gg}

	            assert gg == g

	            g2uFinal[g] = g2u[g].pop()
	            u2gFinal[g2uFinal[g]] = g
	            #tobedelU.add(g2uFinal[g])
	            #tobedelG.add(g)

	        else:
	            #print(u, len(u2g[u]))
	            #solyc id matches uniquely to a uniprot id, but uniprot id matches multiple genes
	            allgenes = list(u2g[u])

	            gstar = allgenes[np.argmax([np.sum(Y[gene2row[gg]]) for gg in allgenes])]

	            for gonidio in allgenes:
	                if gonidio == gstar:
	                    g2uFinal[gonidio] = u
	                    u2gFinal[u] = gonidio
	                else:
	                    tobedelG.add(gonidio)



	    else:
	        nrRev = 0
	        for uu in g2u[g]:
	            if reviewed[uu]:
	                nrRev += 1

	        if nrRev == 0:
	            allgenes = set()
	            for uu in g2u[g]:
	                for gmd in u2g[uu]:
	                    allgenes.add(gmd)

	            if len(allgenes) == 1:
	                gg = allgenes.pop()
	                assert gg == g
	                uu = g2u[g].pop()
	                g2uFinal[g] = uu
	                u2gFinal[uu] = g
	                for gmdu in g2u[g]:
	                    tobedelU.add(gmdu)
	            else:
	                print(g, len(allgenes))
	                for gmd in allgenes:
	                    assert gmd not in g2uFinal
	                #sys.exit(0)

	        elif nrRev == 1:
	            for uu in g2u[g]:
	                if reviewed[uu]:
	                    assert len(u2g[uu]) == 1
	                    g2uFinal[g] = uu
	                    u2gFinal[uu] = g

	                else:
	                    tobedelU.add(uu)

	        else:
	            print(g)

	tobedel = []
	for i, g in enumerate(geneNames):
		if g in g2uFinal or g not in foundHits:
			assert g not in tobedelG
		else:
			tobedel.append(i)

	Y = np.delete(Y, tobedel, axis=0)
	geneNames = np.delete(geneNames, tobedel)

	tobedelTerms = np.where(np.sum(Y,0) == 0)
	Y = np.delete(Y, tobedelTerms, axis=1)
	termNames = np.delete(termNames, tobedelTerms)
	assert np.min(np.sum(Y,1) ) > 0



print(Y.shape)
assert Y.shape == (len(geneNames), len(termNames))

sys.exit(0)

if species == 'tomato':
	with open(path + prefix + 'gene2uniprot.pkl', 'wb') as f:
		pickle.dump(g2uFinal, f)


with open(path + prefix + 'geneNames.pkl', 'wb') as f:
	pickle.dump(geneNames, f)

with open(path + prefix + 'termNames.pkl', 'wb') as f:
	pickle.dump(termNames, f)

save_npz(path + prefix + 'Y.npz', csr_matrix(Y))
