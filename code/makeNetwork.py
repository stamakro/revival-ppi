import pickle
import sys
'''
STEPS
1) map gene id's to uniprot, ignore all other types of ids
2) only use pairs of swissprot entries
3) remove all entries that have isoforms (~90)
4) remove homodimers (protein interacting with itsself)

5) make sure that the sequences of all proteins are known (../sequences/proteins.list)
'''


def makeNetworkYeast():
	path = '../data/yeast/interactions/'

	tair2uni = dict()
	uni2tair = dict()


	with open(path + 'gene2uniprot.tab') as f:
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

	for k in tair2uni:
		assert len(tair2uni[k]) == 1
		assert tair2uni[k][0] in uni2tair

	for k in uni2tair:
		assert len(uni2tair[k]) == 1
		assert uni2tair[k][0] in tair2uni


	for k in tair2uni:
		tair2uni[k] = tair2uni[k][0]


	for k in uni2tair:
		uni2tair[k] = uni2tair[k][0]

	partners = dict()
	with open(path + 'ppi-clean', 'w') as fw:
		with open(path + 'ppi-biogrid-physical.txt') as fr:
			for i, line in enumerate(fr):

				fields = line.split('\t')
				p1 = fields[0]
				p2 = fields[1]


				if p1 == p2:
					#homodimer etc.
					continue


				if p1 not in tair2uni or p2 not in tair2uni:
					continue

				p1 = tair2uni[p1]
				p2 = tair2uni[p2]

				if p1 in partners and p2 in partners[p1]:
					#duplicate line
					#print('duplicate line')
					continue

				if p2 in partners and p1 in partners[p2]:
					#print('reverse')
					#we've seen the reverse order pair
					continue

				if p1 not in partners:
					partners[p1] = set([p2])
				else:
					partners[p1].add(p2)


				fw.write(p1 + '\t' + p2 + '\n')



def makeNetworkArabidopsis():
	path = '../data/arabidopsis/interactions/'

	tair2uni = dict()
	uni2tair = dict()

	with open(path + 'tair2uniprot.tab') as f:
		for i, line in enumerate(f):
			fields = line.split('\t')

			if fields[2] != 'reviewed':
				continue

			uniprot = fields[0]

			geneNames = fields[-1].split()
			for g in geneNames:
				if g.upper()[:4] in ['AT1G', 'AT2G', 'AT3G', 'AT4G', 'AT5G']:
					break
			else:
				if i > 0:
					print(uniprot + ' not there' + str(i))

			gene = g.upper()
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

		for k in tair2uni:
			assert len(tair2uni[k]) == 1
			assert tair2uni[k][0] in uni2tair

		for k in uni2tair:
			assert len(uni2tair[k]) == 1
			assert uni2tair[k][0] in tair2uni


		for k in tair2uni:
			tair2uni[k] = tair2uni[k][0]


		for k in uni2tair:
			uni2tair[k] = uni2tair[k][0]

		partners = dict()
		with open(path + 'ppi-clean', 'w') as fw:
			with open(path + 'ppi-biogrid-physical.txt') as fr:
				for i, line in enumerate(fr):

					fields = line.split('\t')
					p1 = fields[0]
					p2 = fields[1]


					if p1 == p2:
						#homodimer etc.
						continue


					if p1 not in tair2uni or p2 not in tair2uni:
						continue

					p1 = tair2uni[p1]
					p2 = tair2uni[p2]

					if p1 in partners and p2 in partners[p1]:
						#duplicate line
						#print('duplicate line')
						continue

					if p2 in partners and p1 in partners[p2]:
						#print('reverse')
						#we've seen the reverse order pair
						continue

					if p1 not in partners:
						partners[p1] = set([p2])
					else:
						partners[p1].add(p2)


					fw.write(p1 + '\t' + p2 + '\n')


def makeNetworkTomato():

	with open('../data/tomato/annotations/P_geneNames.pkl', 'rb') as f:
		geneNames = pickle.load(f)

	path = '../data/tomato/interactions/'

	etg2solyc = dict()
	solyc2etg = dict()

	with open(path + 'tomato.entrez_2_string.2018.tsv') as f:
		for line in f:
			if line[0] == '#':
				continue

			fields = line.split()

			assert fields[0] == '4081'

			etg = fields[1]
			sol = fields[2].split('.')[1]

			#print etg, sol

			etg = 'ETG' + etg

			if etg in etg2solyc:
				if etg2solyc[etg] in geneNames and sol in geneNames:
					print ('fuck')

				if etg2solyc[etg] in geneNames and sol not in geneNames:
					print(etg, sol)
					continue

			etg2solyc[etg] = sol
			assert sol not in solyc2etg
			solyc2etg[sol] = etg


	etg2solyc['ETG544038'] = 'Solyc06g074350'
	etg2solyc['ETG100736524'] = 'Solyc08g048390'
	etg2solyc['ETG100736519'] = 'Solyc11g020670'
	etg2solyc['ETG100736515'] = 'Solyc02g077250'
	etg2solyc['ETG100736461'] = 'Solyc04006980'
	etg2solyc['ETG100736455'] = 'Solyc06g069240'
	etg2solyc['ETG100736456'] = 'Solyc07g053410'
	etg2solyc['ETG100736247'] = 'Solyc01g103780'
	etg2solyc['ETG543565'] = 'Solyc04g012120'
	partners = dict()


	with open(path + 'ppi-clean', 'w') as fw:
		with open(path + 'ppi-biogrid-physical.txt') as fr:
			for i, line in enumerate(fr):

				fields = line.split('\t')

				p1 = fields[0]
				p2 = fields[1]

				if p1 == p2:
					#homodimer etc.

					continue

				if p1 not in etg2solyc and p1.split('.')[0] != 'Solyc01g094320':

					continue

				if p2 not in etg2solyc and p2.split('.')[0] != 'Solyc01g094320':

					continue

				if p1 in etg2solyc:
					p1 = etg2solyc[p1]
				else:
					p1 = p1.split('.')[0]

				if p2 in etg2solyc:
					p2 = etg2solyc[p2]
				else:
					p2 = p2.split('.')[0]

				assert p1 != p2

				if p1 in partners and p2 in partners[p1]:
					#duplicate line
					continue

				if p2 in partners and p1 in partners[p2]:
					#we've seen the reverse order pair
					continue

				if p1 not in partners:
					partners[p1] = set(p2)
				else:
					partners[p1].add(p2)


				fw.write(p1 + '\t' + p2 + '\n')






species = sys.argv[1]
if species == 'yeast':
	makeNetworkYeast()
elif species == 'arabidopsis':
	makeNetworkArabidopsis()
elif species == 'tomato':
	makeNetworkTomato()
