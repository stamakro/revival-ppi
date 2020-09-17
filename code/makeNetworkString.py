import pickle
import sys
import numpy as np

def makeNetworkYeast(annotationPrefix, excludeUnannotated=True):
	path = '../data/' + species + '/interactions/'


	with open('../data/yeast/annotations/' + annotationPrefix +'geneNames.pkl', 'rb') as f:
		geneNames = pickle.load(f)
	annotatedGenes = set(geneNames)

	tair2uni = dict()
	uni2tair = dict()

	with open(path + 'gene2uniprotString.tab') as f:
	    for i, line in enumerate(f):
	        if i == 0:
	            continue

	        fields = line.split('\t')

	        if fields[1] != 'reviewed':
	            continue

	        uniprot = fields[0]

	        geneNames = fields[-1].split(',')
	        for g in geneNames:

	            gene = g.upper().strip()
	            if gene in tair2uni:
	                tair2uni[gene].append(uniprot)

	            else:
	                tair2uni[gene] = [uniprot]


	            if uniprot in uni2tair:
	                uni2tair[uniprot].append(gene)
	            else:
	                uni2tair[uniprot] = [gene]

	tobedelUni = set()
	tobedelTair = set()

	for k in tair2uni:
	    if len(tair2uni[k]) > 1:
	        tobedelTair.add(k)
	        for uu in tair2uni[k]:
	            tobedelUni.add(uu)

	for i in tobedelTair:
	    del tair2uni[i]

	for i in tobedelUni:
	    del uni2tair[i]

	tobedelUni = set()
	tobedelTair = set()

	for k in uni2tair:
	    if len(uni2tair[k]) > 1:
	        tobedelUni.add(k)
	        for uu in uni2tair[k]:
	            tobedelTair.add(uu)

	for i in tobedelTair:
	    del tair2uni[i]

	for i in tobedelUni:
	    del uni2tair[i]


	tobedelTair = set()
	for k in tair2uni:
	    assert len(tair2uni[k]) == 1
	    if tair2uni[k][0] not in uni2tair:
	        tobedelTair.add(k)

	for kk in tobedelTair:
	    del tair2uni[kk]


	tobedelUni = set()
	for k in uni2tair:
	    assert len(uni2tair[k]) == 1
	    if uni2tair[k][0] not in tair2uni:
	        tobedelUni.add(k)

	for kk in tobedelUni:
	    del uni2tair[kk]

	for k in tair2uni:
	    tair2uni[k] = tair2uni[k][0]


	for k in uni2tair:
	    uni2tair[k] = uni2tair[k][0]


	partners = dict()

	proteins = set()

	with open(path + 'ppi-string.txt') as fr:
	    for i, line in enumerate(fr):
	        if i == 0:
	            fields = line.split()

	            fileNames = fields[2:-1]

	            files = []
	            for ff in fileNames:
	                if excludeUnannotated:
	                    filepointer = open(path + 'features/' + str(ff), 'w')
	                else:
	                    filepointer = open(path + 'features/' + str(ff) + '+unannotated', 'w')
	                files.append(filepointer)


	        else:
	            fields = line.split(' ')
	            p1 = fields[0].split('.')[1]
	            p2 = fields[1].split('.')[1]

	            if p1 == p2:
	                #homodimer etc.
	                continue

	            if p1 not in tair2uni or p2 not in tair2uni:
	                continue


	            p1 = tair2uni[p1]
	            p2 = tair2uni[p2]

	            if (p1 not in annotatedGenes or p2 not in annotatedGenes) and excludeUnannotated:
	                continue

	            assert p1 != p2

	            proteins.add(p1)
	            proteins.add(p2)

	            probs = np.array(fields[2:-1]).astype(int) / 1000.
	            for j, pp in enumerate(probs):
	                if pp > 0:
	                    files[j].write(p1 + '\t' + p2 + '\t' + str(pp) + '\n')

	for ff in files:
	        ff.close()


def makeNetworkArabidopsis(annotationPrefix, excludeUnannotated=True):
	path = '../data/arabidopsis/interactions/'

	with open('../data/arabidopsis/annotations/' + annotationPrefix +'geneNames.pkl', 'rb') as f:
		geneNames = pickle.load(f)
	annotatedGenes = set(geneNames)

	tair2uni = dict()
	uni2tair = dict()

	with open(path + 'gene2uniprotString.tab') as f:
	    for i, line in enumerate(f):
	        if i == 0:
	            continue

	        fields = line.split('\t')

	        if fields[1] != 'reviewed':
	            continue

	        uniprot = fields[0]

	        geneNames = fields[-1].split(' ')

	        for g in geneNames:

	            gene = g.upper().strip()
	            if len(gene) > 4 and gene[:2] == 'AT' and gene[3] == 'G':
	                if gene in tair2uni:
	                	tair2uni[gene].append(uniprot)

	                else:
	                	tair2uni[gene] = [uniprot]

	                if uniprot in uni2tair:
	                	uni2tair[uniprot].append(gene)
	                else:
	                	uni2tair[uniprot] = [gene]


	tobedelUni = set()
	tobedelTair = set()

	for k in tair2uni:
	    if len(tair2uni[k]) > 1:
	        tobedelTair.add(k)
	        for uu in tair2uni[k]:
	            tobedelUni.add(uu)

	for i in tobedelTair:
	    del tair2uni[i]

	for i in tobedelUni:
	    del uni2tair[i]

	tobedelUni = set()
	tobedelTair = set()

	for k in uni2tair:
	    if len(uni2tair[k]) > 1:
	        tobedelUni.add(k)
	        for uu in uni2tair[k]:
	            tobedelTair.add(uu)

	for i in tobedelTair:
	    del tair2uni[i]

	for i in tobedelUni:
	    del uni2tair[i]


	tobedelTair = set()
	for k in tair2uni:
	    assert len(tair2uni[k]) == 1
	    if tair2uni[k][0] not in uni2tair:
	        tobedelTair.add(k)

	for kk in tobedelTair:
	    del tair2uni[kk]


	tobedelUni = set()
	for k in uni2tair:
	    assert len(uni2tair[k]) == 1
	    if uni2tair[k][0] not in tair2uni:
	        tobedelUni.add(k)

	for kk in tobedelUni:
	    del uni2tair[kk]

	for k in tair2uni:
	    tair2uni[k] = tair2uni[k][0]


	for k in uni2tair:
	    uni2tair[k] = uni2tair[k][0]


	partners = dict()

	proteins = set()

	with open(path + 'ppi-string.txt') as fr:
	    for i, line in enumerate(fr):
	        if i == 0:
	            fields = line.split()

	            fileNames = fields[2:-1]

	            files = []
	            for ff in fileNames:
	                if excludeUnannotated:
	                    filepointer = open(path + 'features/' + str(ff), 'w')
	                else:
	                    filepointer = open(path + 'features/' + str(ff) + '+unannotated', 'w')
	                    
	                files.append(filepointer)


	        else:
	            fields = line.split(' ')
	            p1 = fields[0].split('.')[1]
	            p2 = fields[1].split('.')[1]

	            if p1 == p2:
	                #homodimer etc.
	                continue


	            if p1 not in tair2uni or p2 not in tair2uni:
	                continue


	            p1 = tair2uni[p1]
	            p2 = tair2uni[p2]

	            if (p1 not in annotatedGenes or p2 not in annotatedGenes) and excludeUnannotated:
	                continue


	            assert p1 != p2

	            proteins.add(p1)
	            proteins.add(p2)

	            probs = np.array(fields[2:-1]).astype(int) / 1000.
	            for j, pp in enumerate(probs):
	                if pp > 0:
	                    files[j].write(p1 + '\t' + p2 + '\t' + str(pp) + '\n')

	for ff in files:
	        ff.close()



def makeNetworkTomato(annotationPrefix, excludeUnannotated=True):
	path = '../data/' + species + '/interactions/'

	with open('../data/tomato/annotations/' + annotationPrefix +'geneNames.pkl', 'rb') as f:
		geneNames = pickle.load(f)
	annotatedGenes = set(geneNames)


	proteins = set()

	with open(path + 'ppi-string.txt') as fr:
	    for i, line in enumerate(fr):
	        if i == 0:
	            fields = line.split()

	            fileNames = fields[2:-1]

	            files = []
	            for ff in fileNames:
	                if excludeUnannotated:
	                    filepointer = open(path + 'features/' + str(ff), 'w')
	                else:
	                    filepointer = open(path + 'features/' + str(ff) + '+unannotated', 'w')

	                files.append(filepointer)


	        else:
	            fields = line.split(' ')
	            p1 = fields[0].split('.')[1]
	            p2 = fields[1].split('.')[1]

	            if p1 == p2:
	                #homodimer etc.
	                continue

	            if (p1 not in annotatedGenes or p2 not in annotatedGenes) and excludeUnannotated:
	                continue


	            proteins.add(p1)
	            proteins.add(p2)

	            probs = np.array(fields[2:-1]).astype(int) / 1000.
	            for j, pp in enumerate(probs):
	                if pp > 0:
	                    files[j].write(p1 + '\t' + p2 + '\t' + str(pp) + '\n')

	for ff in files:
	        ff.close()



def makeNetworkEcoli(annotationPrefix, excludeUnannotated=True):
	path = '../data/' + species + '/interactions/'


	with open('../data/ecoli/annotations/' + annotationPrefix +'geneNames.pkl', 'rb') as f:
		geneNames = pickle.load(f)
	annotatedGenes = set(geneNames)

	tair2uni = dict()
	uni2tair = dict()

	with open(path + 'ecoli_string2uniprot.tsv') as f:
	    for i, line in enumerate(f):

	        fields = line.split('\t')

	        uni = fields[1].split('|')[0]
	        tair = fields[2].split('.')[1]

	        assert uni not in uni2tair
	        assert tair not in tair2uni

	        uni2tair[uni] = tair
	        tair2uni[tair] = uni

	partners = dict()

	proteins = set()

	with open(path + 'ppi-string.txt') as fr:
	    for i, line in enumerate(fr):
	        if i == 0:
	            fields = line.split()

	            fileNames = fields[2:-1]

	            files = []
	            for ff in fileNames:
	                if excludeUnannotated:
	                    filepointer = open(path + 'features/' + str(ff), 'w')
	                else:
	                    filepointer = open(path + 'features/' + str(ff) + '+unannotated', 'w')

	                files.append(filepointer)


	        else:
	            fields = line.split(' ')
	            p1 = fields[0].split('.')[1]
	            p2 = fields[1].split('.')[1]

	            if p1 == p2:
	                #homodimer etc.
	                continue

	            if p1 not in tair2uni or p2 not in tair2uni:
	                continue


	            p1 = tair2uni[p1]
	            p2 = tair2uni[p2]


	            if (p1 not in annotatedGenes or p2 not in annotatedGenes) and excludeUnannotated:
	                continue


	            assert p1 != p2

	            proteins.add(p1)
	            proteins.add(p2)

	            probs = np.array(fields[2:-1]).astype(int) / 1000.
	            for j, pp in enumerate(probs):
	                if pp > 0:
	                    files[j].write(p1 + '\t' + p2 + '\t' + str(pp) + '\n')

	for ff in files:
	        ff.close()






species = sys.argv[1]

try:
	annotationPrefix = sys.argv[2]
except IndexError:
	annotationPrefix = 'P_'

try:
	excludeUnannotated = bool(int(sys.argv[3]))
except IndexError:
	excludeUnannotated = True

if species == 'yeast':
	makeNetworkYeast(annotationPrefix, excludeUnannotated)
elif species == 'arabidopsis':
	makeNetworkArabidopsis(annotationPrefix, excludeUnannotated)
elif species == 'tomato':
	makeNetworkTomato(annotationPrefix, excludeUnannotated)
elif species == 'ecoli':
	makeNetworkEcoli(annotationPrefix, excludeUnannotated)
