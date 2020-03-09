import numpy as np
import sys
import pickle
from scipy.sparse import csr_matrix


class Gene:

        def __init__(self, name_):

                self.name = name_

                self.functions = set()

                self.predictedFunctions = set()


        def add_function(self, newF):

                if newF not in self.functions:
                        self.functions.add(newF)



        def getFunctions(self):
                return self.functions


        def getName(self):
                return self.name





class Term:
	def __init__(self):
		self.ID = ''
		self.manyIDs = False
		self.altIDs = []

		self.ontology = None

		self.isObsolete = False

		self.isReplaced = False
		self.replacedBy = None
		self.isLeaf = True
		self.parents = []
		self.children = []
		self.associatedGenes = set()
		self.predictedGenes = set()
		self.IC = 0.0
		self.freq = 0.0


	def getCount(self):
		return len(self.associatedGenes)



def merge_annotation_datasets(dag1, dag2):
	new_dag = dag1

	if len(dag1) != len(dag2):
		print ('Different sizes')
		exit()

	for i in range(len(dag1)):
		if dag1[i].ID != dag2[i].ID:
			print ('Different order')
			exit()

		new_dag[i].associatedGenes = (dag1[i].associatedGenes).union(dag2[i].associatedGenes)


	return new_dag



def computePath2Root(dag, ontology):
	path2root = dict()

	if ontology == 'C':
		root = 5895
	elif ontology == 'F':
		root = 1114
	else:
		root = 6

	computePath2RootAux(dag, root, 0, path2root)

	return path2root



def computePath2RootAux(dag, ind, pos, path2root):

		if ind not in path2root:
			path2root[ind] = pos

		elif path2root[ind] < pos:
			path2root[ind] = pos

		for c in dag[ind].children:
			computePath2RootAux(dag, c, pos+1, path2root)




def calculate_IC(dag, ontology):
	new_dag = dag

	for i,t in enumerate(dag):

		if t.ontology != ontology:
			continue

		if len(t.parents) == 0 or t.getCount() == 0:
			t.IC = 0.0

		else:
			parentSet = dag[t.parents[0]].associatedGenes
			for pInd in t.parents:

				parentSet = parentSet.intersection(dag[pInd].associatedGenes)
				#print '%d %d' % (pInd, len(parentSet))

			if len(parentSet) > 0:
				prob = float(t.getCount()) / float(len(parentSet))

				if prob > 1.0:

					#print 'current term'
					print (t.ID)

					prob = 1.0

				t.IC = (-1.0) * np.log2(prob)

			else:
				#print 'Divide by 0'
				#print t.ID

				t.IC = 0.0

		new_dag[i] = t


	return new_dag



def calculate_IC_all_ontologies(dag):
	new_dag = dag

	for i,t in enumerate(dag):

		#if t.ontology != ontology:
		#	continue

		if len(t.parents) == 0 or t.getCount() == 0:
			t.IC = 0.0

		else:
			parentSet = dag[t.parents[0]].associatedGenes
			for pInd in t.parents:

				parentSet = parentSet.intersection(dag[pInd].associatedGenes)
				#print '%d %d' % (pInd, len(parentSet))

			if len(parentSet) > 0:
				prob = float(t.getCount()) / float(len(parentSet))

				if prob > 1.0:

					print('current term')

					exit()


				t.IC = (-1.0) * np.log2(prob)

			else:
				#print 'Divide by 0'
				#print t.ID

				t.IC = 0.0

		new_dag[i] = t


	return new_dag


def find_all_ancestors(dag, tInd):
    ancestors = set(dag[tInd].parents)

    for pInd in dag[tInd].parents:
        ancestors = ancestors.union(find_all_ancestors(dag, pInd))

    return ancestors


def find_ontology_ancestors(dag, tInd, ontology):
    ancestors = set([ind for ind in dag[tInd].parents if dag[ind].ontology == ontology])
    #ancestors = set(dag[tInd].parents)

    for pInd in dag[tInd].parents:
        ancestors = ancestors.union(find_ontology_ancestors(dag, pInd, ontology))

    return ancestors



def find_common_ancestors_denovo(dag, t1Ind, t2Ind, ontology, all):

    if not all:
        t1ancestors = find_ontology_ancestors(dag, t1Ind, ontology)
        t2ancestors = find_ontology_ancestors(dag, t2Ind, ontology)
    else:
        t1ancestors = find_all_ancestors(dag, t1Ind)
        t2ancestors = find_all_ancestors(dag, t2Ind)


    return t1ancestors.intersection(t2ancestors)


def get_all_frequencies(dag, ontology):

    individual_probability = dict()
    for i, t in enumerate(dag):
        if t.ontology == ontology:
          individual_probability[i] = float(t.getCount()) / float(nrPresentGenes)

    return individual_probability



def dyn_ancestors(dag, ontology, anc):

	leaves = [i for i in range(len(dag)) if dag[i].ontology == ontology and len(dag[i].children) == 0]

	ontoTerms = [i for i in range(len(dag)) if dag[i].ontology == ontology]

	isDone = np.zeros((len(dag),1), int)

	for leaf in leaves:
		recAncestors(dag, anc, leaf, isDone)

	return anc


def recAncestors(dag, anc, i, isDone):

	if isDone[i]:
		return anc[i]


	anc[i] = set(dag[i].parents)



	for pInd in dag[i].parents:
			anc[i] = anc[i].union(recAncestors(dag, anc, pInd, isDone))

	isDone[i] = 1
	#print '!'

	return anc[i]


def getLCA(dag, ontology, ancestors, ti, tj):

	anci = ancestors[ti]

	ancj = ancestors[tj]

	if ti in ancj:
		ms = ti

	elif tj in anci:
		ms = tj

	else:
		common = anci.intersection(ancj)

		if len(common) == 1:
			ms = common.pop()

		else:

			tmpMs = -1.0
			ind = -1

			for t in common:
				if dag[t].IC > tmpMs:
					tmpMs = dag[t].IC
					ind = t

			ms = ind


	return ms




def lin(dag, ontology):
	#dag should be deep-copied at call


	ontoTerms = [i for i in range(len(dag)) if dag[i].ontology == ontology]

	setTerms = set(ontoTerms)

	ancestors = [set() for i in range(len(dag))]

	'''
	if ontology == 'C':
		root = 5895
	elif ontology == 'F':
		root = 1114
	else:
		root = 6
	'''

	dyn_ancestors(dag, ontology, ancestors)


	#return ancestors

	eps = 1e-5
	#test for 0

	lowestCommonAnc = np.zeros((len(ontoTerms), len(ontoTerms)), int)

	linMatrix = np.zeros((len(ontoTerms), len(ontoTerms)), float)

	lll = len(ontoTerms)


	for i in range(len(ontoTerms)):

		#if i % 100 == 0:
		#	print '%d out of %d' % (i, lll)

		lowestCommonAnc[i,i] = i


		pi = dag[ontoTerms[i]].IC

		if pi < eps:
			continue

		linMatrix[i,i] = 1.0

		for j in range(i + 1, len(ontoTerms)):

			pj = dag[ontoTerms[j]].IC

			if pj < eps:
				continue


			lca = getLCA(dag, ontology, ancestors, ontoTerms[i], ontoTerms[j])

			lowestCommonAnc[i, j] = lca

			lsij = 2.0 * dag[lca].IC / (pi + pj)

			linMatrix[i,j] = lsij

			linMatrix[j,i] = lsij


	return [ontoTerms, linMatrix]



def annotate_gene_with_term(gene, dag, termInd):

	if dag[termInd].ID not in gene.getFunctions():

		gene.add_function(dag[termInd].ID)

		dag[termInd].associatedGenes.add(gene)

	for pInd in dag[termInd].parents:

		l =  annotate_gene_with_term(gene, dag, pInd)
		dag = l[0]
		gene = l[1]


	return [dag, gene]

def annotate_gene_with_term2(gene, dag, termInd):

	gene.add_function(dag[termInd].ID)

	dag[termInd].associatedGenes.add(gene)

	for pInd in dag[termInd].parents:

		l =  annotate_gene_with_term2(gene, dag, pInd)
		dag = l[0]
		gene = l[1]


	return [dag, gene]



def read_go_per_gene_ec(fil, obofile=None, excludedEcs=[], allowedExplanations=dict()):

	if len(excludedEcs) == 0:
		assert len(allowedExplanations) == 0

	if obofile is None:
		ontology = read_ontology_from_file('/tudelft.net/staff-bulk/ewi/insy/DBL/smakrod/go/terms.obo')
	else:
		ontology = read_ontology_from_file(obofile)

	tree = ontology[0]
	mapping = ontology[1]


	dupl = [i for i in range(len(tree)) if tree[i].manyIDs]


	allGenes = []

	currentGene = None

	with open(fil) as f:
		for line in f:
			if line[0]== '!':
				continue

			fields = line.split('\t')

			if fields[3] == 'NOT':
				continue


			name = fields[1]
			term = fields [4]

			ec = fields[6]
			ecExplanation = fields[7].split(':')[0]

			skip = False
			if ec  in excludedEcs:
				skip = True
			if ec in allowedExplanations and ecExplanation in allowedExplanations[ec]:
				skip = False

			if skip:
				continue



			stt = term[:3]
			if stt != 'GO:':
				print ('fuck')

			if term not in mapping:
				found = False
				i = 0

				while i < len(dupl) and not found:

					if term in tree[dupl[i]].altIDs:
						term = tree[dupl[i]].ID
						found = True

					i += 1


				if not found:
					print ('still not found: ', end='')
					print (term)


			if currentGene is None:

				currentGene = Gene(name)


				l = annotate_gene_with_term(currentGene, tree, mapping[term])

				tree = l[0]
				currentGene = l[1]



			elif name == currentGene.getName():
				l = annotate_gene_with_term(currentGene, tree, mapping[term])

				tree = l[0]
				currentGene = l[1]


			else:
				allGenes.append(currentGene)

				currentGene = Gene(name)


				l = annotate_gene_with_term(currentGene, tree, mapping[term])

				tree = l[0]
				currentGene = l[1]




	return [allGenes, tree, mapping]



def read_ontology_from_file(fil):

	originalIDs = set()
	obsolete = set()

	mapping = dict()

	useless = set(['\n', 'disjoint_from', 'comment', 'consider', 'created_by', 'creation_date', 'def', 'name', 'property_value', 'subset', 'synonym', 'xref', 'intersection_of'])

	usefullRelations = set(['part_of', 'regulates', 'positively_regulates', 'negatively_regulates', 'occurs_in'])

	dag = []

	with open(fil) as f:
		for line in f:
			#print line
			if line == '\n':
				if dag[currentInd].isObsolete == True:
					obsolete.add(dag[currentInd].ID)
					del dag[currentInd]



			if line.rstrip('\n') == '[Term]':
				t = Term()

				dag.append(t)

				currentInd = len(dag) - 1

			else:
				fields = line.split(':')

				pr = fields[0]


				if pr in useless:
					continue


				elif pr == 'id':

					ID = 'GO:' + fields[2].rstrip('\n')

					originalIDs.add(ID)

					if ID not in mapping:
						'''new term'''
						dag[currentInd].ID = ID
						mapping[t.ID] = currentInd

					else:
						'''seen before as someone's parent or replacement'''
						'''remove newly added object and get right position in the list'''
						del dag[currentInd]
						currentInd = mapping[ID]



				elif pr == 'namespace':
					if fields[1].strip() == 'cellular_component':
						dag[currentInd].ontology = 'C'

					elif fields[1].strip() == 'molecular_function':
						dag[currentInd].ontology = 'F'


					elif fields[1].strip() == 'biological_process':
						dag[currentInd].ontology = 'P'

					else:
						print ('DANGER: unknown ontollogy')


				elif pr == 'alt_id':
					dag[currentInd].manyIDs = True
					newID = 'GO:' + fields[2].rstrip('\n')

					dag[currentInd].altIDs.append(newID)
					mapping[newID] = currentInd


				elif pr == 'is_obsolete':
					if fields[1].strip() == 'true':
						dag[currentInd].isObsolete = True
					else:
						print (fields[1])


				elif pr == 'replaced_by':
					continue
					''' if obsolete, just delete it. obsolete terms do not have parents   '''


				elif pr == 'is_a':
					parentID = 'GO:' + fields[2].split('!')[0].strip()

					if parentID not in mapping:
						tt = Term()
						tt.ID = parentID
						dag.append(tt)

						mapping[parentID] = len(dag) - 1

					if mapping[parentID] not in dag[currentInd].parents:
						'''this is because some are and regulate'''
						dag[currentInd].parents.append( mapping[parentID])
						dag[mapping[parentID]].children.append(currentInd)

					dag[mapping[parentID]].isLeaf = False

				elif pr == 'relationship':

					relType = fields[1].split()[0]

					if relType not in usefullRelations:
						continue

					parentID = 'GO:' + fields[2].split('!')[0].strip()

					if parentID not in mapping:
						tt = Term()
						tt.ID = parentID
						dag.append(tt)

						mapping[parentID] = len(dag) - 1

					if mapping[parentID] not in dag[currentInd].parents:
						'''this is because some are and regulate'''
						dag[currentInd].parents.append( mapping[parentID])
						dag[mapping[parentID]].children.append(currentInd)

					dag[mapping[parentID]].isLeaf = False


				else:
					print (pr)

	return [dag, mapping]

def getLabelMatrix(genes, dag, ontology):

	#this function removes terms with no annotated proteins
	#CAREFUL! After calling it, the dag[x].parents, dag[x].children are no longer valid!


	if not isinstance(next(iter(genes)), str):
		usedGenesNames = [g.name for g in genes]
	else:
		usedGenesNames = genes


	[usedTermsIndices, usedTermsNames] = getUsedTerms(dag, ontology, set(usedGenesNames))


	geneDict = dict()

	for i, t in enumerate(usedGenesNames):
		geneDict[t] = i

	data = np.zeros((len(genes), len(usedTermsIndices)), int)

	for j, tInd in enumerate(usedTermsIndices):
		assocGenes = [g.name for g in dag[tInd].associatedGenes]

		for gg in assocGenes:
			try:
				data[geneDict[gg], j] = 1

			except KeyError:
				continue


	return [data, usedTermsNames, usedGenesNames]



def getUsedTerms(dag, ontology, targetGenes):

	used = []
	names = []
	for i, t in enumerate(dag):

		if t.ontology == ontology:

			assocG = set([g.name for g in t.associatedGenes])


			if len(assocG.intersection(targetGenes)) > 0:
				used.append(i)
				names.append(t.ID)
	return [used, names]



def getMaxPathToRoot(termNames, ontology, dag=None):

	if dag is None:
		dag, mapping = fr.read_ontology_from_file('/tudelft.net/staff-bulk/ewi/insy/DBL/smakrod/go/terms.obo')

	allpaths = fr.computePath2Root(dag, ontology)

	pathLengths = [allpaths[mapping[t]] for t in termNames]

	return np.array(pathLengths)


def getParentsCoord(termNames, ontology, dag=None, mapping=None):
    if dag is None:
        dag, mapping = fr.read_ontology_from_file('/tudelft.net/staff-bulk/ewi/insy/DBL/smakrod/go/terms.obo')
    if ontology == 'P':
        root = 'GO:0008150'
    else:
        if ontology == 'F':
            root = 'GO:0003674'
        else:
            print ('wrong ontology')
            sys.exit(1)
    parentsCoord = dict()
    for i, tname in enumerate(termNames):
        #print tname
        dagInd = mapping[tname]
        parentsCoord[i] = []
        for pInd in dag[dagInd].parents:
            if dag[pInd].ontology != ontology or dag[pInd].ID == root:
                continue
            parLoc = termNames.index(dag[pInd].ID)
            parentsCoord[i].append(parLoc)

    return parentsCoord


def removeEmptyProteins(X, Y):
    nrTerms = np.sum(Y, 1)
    zeros = np.where(nrTerms == 0)
    X = np.delete(X, zeros, axis=0)
    X = np.delete(X, zeros, axis=1)
    Y = np.delete(Y, zeros, axis=0)
    return [
     X, Y, zeros[0]]


def removeEmptyTestProteins(X, Y):
    nrTerms = np.sum(Y, 1)
    zeros = np.where(nrTerms == 0)
    X = np.delete(X, zeros, axis=0)
    Y = np.delete(Y, zeros, axis=0)
    return [
     X, Y]


def removeRareTermsOnly(Y, G, label_fraction):
    label_fraction /= 100.0
    nrProteins, nrTerms = Y.shape
    fractions = np.sum(Y, 0) / float(nrProteins)
    tokeep = np.where(fractions > label_fraction)[0]
    Y = Y[:, tokeep]
    G = G[:, tokeep][tokeep, :]
    return [
     Y, G]


def removeRareTerms(X, Y, label_fraction):
    label_fraction /= 100.0
    nrProteins, nrTerms = Y.shape
    fractions = np.sum(Y, 0) / float(nrProteins)
    tokeep = np.where(fractions > label_fraction)[0]
    Y = Y[:, tokeep]
    X, Y, removedProteins = removeEmptyProteins(X, Y)
    return [
     X, Y, tokeep, removedProteins]


def calculateIC(Ytrain, parentsCoord):


    #root = np.where(np.sum(Ytrain,0) == Ytrain.shape[0])[0]

    #if len(root) != 1:
    #    return calculateICnoRoot(Ytrain, parentsCoord)

    nrProteins, nrTerms = Ytrain.shape
    ic = np.zeros((nrTerms,), float)
    for i in range(nrTerms):
        numerator = Ytrain[:, i]
        denominator = np.ones(nrProteins)
        for p in parentsCoord[i]:
            denominator = np.minimum(denominator, Ytrain[:, p])

        if np.sum(denominator) > 0 and np.sum(numerator) > 0:
            p = float(np.sum(numerator)) / float(np.sum(denominator))
            if p < 0.0 or p > 1.0:
                print ('wtf')
                print (i, p)
            ic[i] = -1.0 * np.log2(p)
        else:
            ic[i] = 0.0

    return ic



def normalizedRemainingUncertainty(Ytrue, Ypred, termIC, avg=False):
    num =  np.logical_and(Ytrue == 1, Ypred == 0).astype(float).dot(termIC)
    denom =  np.logical_or(Ytrue == 1, Ypred == 1).astype(float).dot(termIC)

    nru = num / denom

    if avg:
        nru = np.mean(nru)

    return nru


def normalizedMisInformation(Ytrue, Ypred, termIC, avg=False):
    num =  np.logical_and(Ytrue == 0, Ypred == 1).astype(float).dot(termIC)
    denom =  np.logical_or(Ytrue == 1, Ypred == 1).astype(float).dot(termIC)

    nmi = num / denom

    if avg:
        nmi = np.mean(nmi)

    return nmi


def normalizedSemanticDistance(Ytrue, Ypred, termIC, avg=True):
    ru = normalizedRemainingUncertainty(Ytrue, Ypred, termIC, avg)
    mi = normalizedMisInformation(Ytrue, Ypred, termIC, avg)
    sd = np.sqrt(ru ** 2 + mi ** 2)

    #if avg:
    #    ru = np.mean(ru)
    #    mi = np.mean(mi)
    #    sd = np.mean(sd)

    return [ru, mi, sd]


def remainingUncertainty(Ytrue, Ypred, termIC, avg=False):
    ru = np.logical_and(Ytrue == 1, Ypred == 0).astype(float).dot(termIC)
    if avg:
        ru = np.mean(ru)

    return ru
    '''
    nrProteins, nrTerms = Ytrue.shape
    
    ru = np.zeros((nrProteins,), float)
    for i in range(nrProteins):
        pr = Ypred[i]
        gt = Ytrue[i]
        fn = np.intersect1d(np.where(gt == 1), np.where(pr == 0))
        for fni in fn:
            ru[i] += termIC[fni]

    if avg:
        ru = np.mean(ru)
    return ru
    '''

def misInformation(Ytrue, Ypred, termIC, avg=False):
    mi = np.logical_and(Ytrue == 0, Ypred == 1).astype(float).dot(termIC)
    '''
    nrProteins, nrTerms = Ytrue.shape
    mi = np.zeros((nrProteins,), float)
    for i in range(nrProteins):
        pr = Ypred[i]
        gt = Ytrue[i]
        fp = np.intersect1d(np.where(gt == 0), np.where(pr == 1))
        for fpi in fp:
            mi[i] += termIC[fpi]
    '''
    if avg:
        mi = np.mean(mi)
    return mi


def semanticDistance(Ytrue, Ypred, termIC):
    ru = remainingUncertainty(Ytrue, Ypred, termIC, True)
    mi = misInformation(Ytrue, Ypred, termIC, True)
    sd = np.sqrt(ru ** 2 + mi ** 2)
    return [ru, mi, sd]


def transform(X):
    m = float(np.amin(X))
    M = float(np.amax(X))
    return (X - m) / (M - m)


def calculateGOseq(ontology):
    organism = 'thalia'
    dag, mapping = fr.read_ontology_from_file('/tudelft.net/staff-bulk/ewi/insy/DBL/smakrod/go/terms.obo')
    Y, termNames, _ = fr.getDataLabels(organism, True, ontology, 8677)
    go = np.eye(len(termNames))
    for i, tname in enumerate(termNames):
        dagInd = mapping[tname]
        for pInd in dag[dagInd].parents:
            if dag[pInd].ontology != ontology:
                continue
            parLoc = termNames.index(dag[pInd].ID)
            go[(i, parLoc)] = 1.0
            go[(parLoc, i)] = 1.0

    print (go.shape)
    go = csr_matrix(go)
    with open('go' + ontology + '.pkl', 'wb') as (f):
        pickle.dump(go, f)
    return go


def calculateGOallseq(ontology):
    organism = 'thalia'
    dag, mapping = fr.read_ontology_from_file('/tudelft.net/staff-bulk/ewi/insy/DBL/smakrod/go/terms.obo')
    Y, termNames, _ = fr.getDataLabels(organism, True, ontology, 8000)
    go = DFS(dag, termNames, mapping, ontology)
    go = csr_matrix(go)
    with open('go_allancestors' + ontology + '.pkl', 'wb') as (f):
        pickle.dump(go, f)
    return go


def DFS(dag, termNames, mapping, ontology):
    if ontology == 'C':
        root = 5895
    else:
        if ontology == 'F':
            root = 1114
        else:
            root = 6
    gomat = np.eye(len(termNames))
    gomat2, uselesspath = visitNode(gomat, dag, termNames, mapping, [], termNames.index(dag[root].ID))
    return gomat2


def visitNode(gomat, dag, termNames, mapping, path, node):
    path.append(node)
    for cc in path:
        gomat[(node, cc)] = 1
        gomat[(cc, node)] = 1

    for ch in dag[mapping[termNames[node]]].children:
        if dag[ch].ID not in termNames:
            continue
        gomat, path = visitNode(gomat, dag, termNames, mapping, path, termNames.index(dag[ch].ID))

    del path[-1]
    return [
     gomat, path]


def getLCAs(termNames, ontology, dag=None, mapping=None):
    if dag is None:
        dag, mapping = fr.read_ontology_from_file('/tudelft.net/staff-bulk/ewi/insy/DBL/smakrod/go/terms.obo')
    if ontology == 'P':
        root = 'GO:0008150'
    else:
        if ontology == 'F':
            root = 'GO:0003674'
        else:
            print ('wrong ontology')
            sys.exit(1)
    ancestors = dict()
    for i, t in enumerate(termNames):
        ancestors[t] = set()
        for pInd in dag[mapping[t]].parents:
            getAncestorsRecursively(dag, pInd, ancestors[t], ontology)

        ancestors[t].remove(root)

    lcas = dict()
    for i, t1 in enumerate(termNames):
        print (i)
        for j, t2 in enumerate(termNames):
            if i == j:
                lcas[(t1, t2)] = t1

            elif t1 in ancestors[t2]:
                lcas[(t1, t2)] = t1

            elif t2 in ancestors[t1]:
                lcas[(t1, t2)] = t2

            else:
                cas = ancestors[t1].intersection(ancestors[t2])
                for t in cas:
                    others = cas.difference(set([t]))
                    if others == ancestors[t]:
                        lcas[(t1, t2)] = t
                        break

    print ('save')
    with open('lowestCommonAncestors.pkl', 'wb') as (f):
        pickle.dump(lcas, f)
    with open('ancestors.pkl', 'wb') as (f):
        pickle.dump(ancestors, f)
    return


def getAncestorsRecursively(dag, current, ancestors, ontology):
    if dag[current].ontology == ontology:
        ancestors.add(dag[current].ID)
    for pInd in dag[current].parents:
        getAncestorsRecursively(dag, pInd, ancestors, ontology)


def getLCAsCAFA(termNames, ontology, dag=None, mapping=None):
    if dag is None:
        dag, mapping = fr.read_ontology_from_file('/tudelft.net/staff-bulk/ewi/insy/DBL/smakrod/go/terms_09_17.obo')
    if ontology == 'P':
        root = 'GO:0008150'
    else:
        if ontology == 'F':
            root = 'GO:0003674'
        else:
            print ('wrong ontology')
            sys.exit(1)
    ancestors = dict()
    for i, t in enumerate(termNames):
        ancestors[t] = set()
        for pInd in dag[mapping[t]].parents:
            getAncestorsRecursively(dag, pInd, ancestors[t], ontology)

        ancestors[t].remove(root)

    lcas = dict()
    for i, t1 in enumerate(termNames):
        print (i)
        for j, t2 in enumerate(termNames):
            if i == j:
                lcas[(t1, t2)] = t1
            else:
                cas = ancestors[t1].intersection(ancestors[t2])
                for t in cas:
                    others = cas.difference(set([t]))
                    if others == ancestors[t]:
                        lcas[(t1, t2)] = t
                    break

    with open('lowestCommonAncestors.pkl', 'wb') as (f):
        pickle.dump(lcas, f)
    with open('ancestors.pkl', 'wb') as (f):
        pickle.dump(ancestors, f)
    return


def calculateGOseqCAFA(ontology, termNames):
    organism = 'thalia'
    dag, mapping = fr.read_ontology_from_file('/tudelft.net/staff-bulk/ewi/insy/DBL/smakrod/go/terms_09_17.obo')
    if ontology == 'P':
        root = 'GO:0008150'
    else:
        if ontology == 'F':
            root = 'GO:0003674'
        else:
            print ('wrong ontology')
            sys.exit(1)
    go = np.eye(len(termNames))
    for i, tname in enumerate(termNames):
        dagInd = mapping[tname]
        for pInd in dag[dagInd].parents:
            if dag[pInd].ontology != ontology or dag[pInd] == root:
                continue
            parLoc = termNames.index(dag[pInd].ID)
            go[(i, parLoc)] = 1.0
            go[(parLoc, i)] = 1.0

    go = csr_matrix(go)
    with open('go' + ontology + '.pkl', 'wb') as (f):
        pickle.dump(go, f)
    return go
