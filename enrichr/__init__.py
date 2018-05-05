VERSION = '0.0.0dev'

import requests, json, re
from os import path
from collections import namedtuple

# types
Enrichr_Listid   = namedtuple('Enrichr_Listid', ['userListId', 'shortId'])
Enrichr_Genelist = namedtuple('Enrichr_Genelist', ['genes', 'description'])
Enrichr_Term     = namedtuple('Enrichr_Term', [
	'rank',
	'term',
	'pval',
	'zscore',
	'combinedScore',
	'overlappingGenes',
	'adjPval',
	'oldPval',
	'oldAdjPval'
])
Enrichr_Library = namedtuple('Enrichr_Library', [
	'name',
	'category',
	'hasGrid',
	'isFuzzy',
	'format',
	'description',
	'terms'
])

class Enrichr(object):
	# urls
	URL_ADDLIST = 'http://amp.pharm.mssm.edu/Enrichr/addList'
	URL_VIEW    = 'http://amp.pharm.mssm.edu/Enrichr/view?userListId={listid}'
	URL_ENRICH  = 'http://amp.pharm.mssm.edu/Enrichr/enrich?userListId={listid}&backgroundType={library}'
	URL_GENEMAP = 'http://amp.pharm.mssm.edu/Enrichr/genemap?json=true&setup=true&gene={gene}'
	URL_EXPORT  = 'http://amp.pharm.mssm.edu/Enrichr/export'

	# msgs
	MSG_NO_LISTID          = 'You have to successfully add a list first'
	MSG_NO_LIBRARY         = 'You have to set a library'
	MSG_NO_RESULTS         = 'You have to run an enrich first'
	MSG_ERROR_ADD_GENELIST = 'Error adding gene list'
	MSG_ERROR_GET_GENELIST = 'Error getting gene list'
	MSG_ERROR_ENRICH       = 'Error fetching enrichment results'
	MSG_ERROR_SEARCH_TERM  = 'Error searching for terms'

	def __init__(self, library = 'KEGG_2016'):
		self.listid  = None
		self.library = library
		self.results = []

	def addListFromFile(self, genefile, col = 0, delimit = '\t', skip = 0):
		description  = path.basename(genefile)
		genes = set()
		with open(genefile) as f:
			for i in range(skip):
				f.readline()
			for line in f:
				line = line.strip()
				if not line: continue
				genes.add(line.split('\t')[col])
		return self.addList(genes, description)

	def addList(self, genes, description = 'Enrichr gene list.'):
		genes_str = '\n'.join(set(genes))
		payload = {
			'list': (None, genes_str),
			'description': (None, description)
		}
		response = requests.post(Enrichr.URL_ADDLIST, files=payload)
		if not response.ok:
			msg = re.sub('<[^<]+?>', '', response.text).splitlines()
			#msg = response.text.splitlines()
			msg = [line.strip() for line in msg if line.strip()]
			msg = ', '.join(msg)
			raise Exception(Enrichr.MSG_ERROR_ADD_GENELIST + ': ' + msg)

		data = json.loads(response.text)
		self.listid = data['userListId']
		return Enrichr_Listid(**data)

	def view(self, listid = None):
		listid = listid or self.listid
		if not listid:
			raise Exception(Enrichr.MSG_NO_LISTID)
		self.listid = listid

		response = requests.get(Enrichr.URL_VIEW.format(listid = listid))
		if not response.ok:
			raise Exception(Enrichr.MSG_ERROR_GET_GENELIST)

		data = json.loads(response.text)
		return Enrichr_Genelist(genes = data['genes'], description = data['description'])

	def enrich(self, library = None, listid = None):
		listid = listid or self.listid
		if not listid:
			raise Exception(MSG_NO_LISTID)
		self.listid = listid

		library = library or self.library
		if not library:
			raise Exception(MSG_NO_LIBRARY)
		self.library = library

		response = requests.get(Enrichr.URL_ENRICH.format(listid = listid, library = library))
		if not response.ok:
			raise Exception(Enrichr.MSG_ERROR_ENRICH)
		data = json.loads(response.text)
		self.results = sorted([Enrichr_Term._make(dt) for dt in data[library]], key = lambda term: float(term.adjPval))
		return self.results

	def genemap(self, gene):
		response = requests.get(Enrichr.URL_GENEMAP.format(gene = gene))
		if not response.ok:
			raise Exception(Enrichr.MSG_ERROR_SEARCH_TERM)

		data = json.loads(response.text)
		descriptions = {desc['name']:(desc['description'] if 'description' in desc else '') for desc in data['descriptions']}
		for category in data['categories']:
			for lib in category['libraries']:
				yield Enrichr_Library(
					name        = lib['name'],
					category    = category['name'],
					hasGrid     = lib['hasGrid'],
					isFuzzy     = lib['isFuzzy'],
					format      = lib['format'],
					description = descriptions[lib['name']],
					terms       = data['gene'][lib['name']]
				)

	def export(self, outfile, top = None):
		if not self.results:
			raise Exception(Enrichr.MSG_NO_RESULTS)

		results = self.results if top is None else self.results[:top]
		with open(outfile, 'w') as fout:
			fout.write('\t'.join(Enrichr_Term._fields) + '\n')
			for term in results:
				fout.write('\t'.join([
					'|'.join(term.overlappingGenes) if field == 'overlappingGenes' else str(getattr(term, field)) for field in term._fields
				]) + '\n')

	def plot(self, outfile, title = 'Gene enrichment: {library}', top = 10):
		import matplotlib, math
		matplotlib.use('Agg')
		from matplotlib import pyplot as plt
		from matplotlib import gridspec
		from matplotlib import patches

		# already sorted
		results = self.results[:top]

		gs = gridspec.GridSpec(1, 2, width_ratios=[3, 7])
		rownames = [r.term if len(r.term)<40 else r.term[:36] + ' ...' for r in results]
		rnidx    = range(len(rownames))
		ax1 = plt.subplot(gs[0])
		plt.title (title.format(library = self.library), fontweight='bold')

		ax1.xaxis.grid(alpha=.6, ls = '--', zorder = -99)
		plt.subplots_adjust(wspace=.01, left=0.5)
		ax1.barh(rnidx, [len(r.overlappingGenes) for r in results], color='blue', alpha=.6)
		plt.yticks (rnidx, rownames)
		ax1.yaxis.set_ticks_position('none')
		ax1.tick_params(axis='x', colors='blue')
		ax1.spines['top'].set_visible(False)
		ax1.spines['left'].set_visible(False)
		ax1.spines['right'].set_visible(False)
		ax1.spines['bottom'].set_linewidth(1)
		ax1.invert_xaxis()
		ax1.invert_yaxis()
		xticks = ax1.xaxis.get_major_ticks()
		xticks[0].label1.set_visible(False)

		ax2 = plt.subplot(gs[1])
		ax2.xaxis.grid(alpha=.6, ls = '--', zorder = -99)
		ax2.barh(rnidx, [-math.log(r[2], 10) for r in results], color='red', alpha = .6)
		for i, r in enumerate(results):
			t  = str("%.2E" % r[2])
			tx = 0.1
			ty = i + 0.1
			ax2.text(tx, ty, t, fontsize=8)
		ax2.tick_params(axis='x', colors='red')
		ax2.spines['top'].set_visible(False)
		ax2.spines['left'].set_visible(False)
		ax2.spines['right'].set_visible(False)
		ax2.spines['bottom'].set_linewidth(1)
		ax2.invert_yaxis()
		plt.yticks([])
		ng_patch = patches.Patch(color='blue', alpha=.6, label='# overlapping genes')
		pv_patch = patches.Patch(color='red', alpha = .6, label='-log(p-value)')
		plt.figlegend(handles=[ng_patch, pv_patch], labels=['# overlapping genes', '-log(p-value)'], loc="lower center", ncol=2, edgecolor="none")
		plt.savefig(outfile, dpi=300)
