VERSION = '0.0.0dev'

import requests, json, re
from os import path
from subprocess import check_output
from collections import namedtuple
from tempfile import NamedTemporaryFile

# types
Enrichr_Listid   = namedtuple('Enrichr_Listid', ['userListId', 'shortId'])
Enrichr_Genelist = namedtuple('Enrichr_Genelist', ['genes', 'description'])
Enrichr_Term     = namedtuple('Enrichr_Term', [
	'Term',
	'Overlap',
	'Pval',
	'AdjPval',
	'OldPval',
	'OldAdjPval',
	'Z',
	'CombinedScore',
	'Genes'
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
	URL_EXPORT  = 'http://amp.pharm.mssm.edu/Enrichr/export?userListId={listid}&backgroundType={library}'

	# msgs
	MSG_NO_LISTID          = 'You have to successfully add a list first'
	MSG_NO_LIBRARY         = 'You have to set a library'
	MSG_NO_RESULTS         = 'You have to run an enrich first'
	MSG_ERROR_ADD_GENELIST = 'Error adding gene list'
	MSG_ERROR_GET_GENELIST = 'Error getting gene list'
	MSG_ERROR_ENRICH       = 'Error fetching enrichment results'
	MSG_ERROR_SEARCH_TERM  = 'Error searching for terms'

	def __init__(self, library = 'KEGG_2016', cutoff = 0.05, top = 20, Rscript = 'Rscript'):
		"""
		Constructor
		@params:
			`library`: The library to do enrichment against. Default: `KEGG_2016`
				- See all available libraries: http://amp.pharm.mssm.edu/Enrichr/#stats
			`cutoff` : The cutoff for the enriched terms. Default: `0.05`
				- A float indicating cutoff by "AdjPval", or
				- A dict specifying the `value` and `by`
			`top`    : Top N terms. Filtering terms with `cutoff`, whichever comes first. Default: `20`
				- Filters (cutoff and top) are for export and plot
		"""
		# save the session
		self.listid  = None
		self.size    = None
		self.library = library
		self.results = []
		self.cutoff  = cutoff if isinstance(cutoff, dict) else {"by": "AdjPval", "value": cutoff}
		self.top     = top
		self.Rscript = Rscript

	def addListFromFile(self, genefile, col = 0, delimit = '\t', skip = 0):
		"""
		Add a gene list from files
		@params:
			`genefile`: The file containing genes
			`col`     : The column index where the genes located. Default: `0`
			`delimit` : The delimiter for the columns. Default : `\t`
			`skip`    : Skip first N lines. Default : `0`
		@returns:
			A `Enrichr_Listid` with the session ids.
		"""
		description  = path.basename(genefile)
		genes = set()
		with open(genefile) as f:
			for i in range(skip):
				f.readline()
			for line in f:
				line = line.strip()
				if not line: continue
				genes.add(line.split('\t')[col].strip())
		return self.addList(genes, description)

	def addList(self, genes, description = 'Enrichr gene list.'):
		"""
		Add a gene list
		@params:
			`genes`: The genes.
			`description`: The description of the gene list
		@returns:
			A `Enrichr_Listid` with the session ids.
		"""
		self.size = len(genes)
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
		"""
		View added gene list
		@params:
			`listid`: The id of the gene list. Default: `None` (current list)
		@returns:
			`Enrichr_Genelist`
		"""
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
		"""
		Do the enrichment
		@params:
			`library`: THe library to do enrichment against. Default: `None` (current library)
			`listid` : The list id of the session. Default: `None` (current session)
		@returns:
			A list of `Enrichr_Term` (sorted by `adjPval`)
		"""
		listid = listid or self.listid
		if not listid:
			raise Exception(MSG_NO_LISTID)
		self.listid = listid

		library = library or self.library
		if not library:
			raise Exception(MSG_NO_LIBRARY)
		self.library = library

		#response = requests.get(Enrichr.URL_ENRICH.format(listid = listid, library = library))
		#if not response.ok:
		#	raise Exception(Enrichr.MSG_ERROR_ENRICH)
		# doesn't have overlap numbers
		#data = json.loads(response.text)
		#self.results = sorted([Enrichr_Term._make(dt) for dt in data[library]], key = lambda term: float(term.adjPval))
		response = requests.get(Enrichr.URL_EXPORT.format(listid = listid, library = library), stream = True)
		self.results = []
		for line in response.iter_lines():
			if line.startswith('Term'):
				continue
			data = line.split('\t')
			data[2:-1] = [float(d) for d in data[2:-1]]
			data[0] = data[0].replace('_Homo sapiens', '')
			self.results.append(Enrichr_Term._make(data))
		self.results = sorted(self.results, key = lambda x: x.AdjPval)
		return self.results

	def genemap(self, gene):
		"""
		Find terms that contain a given gene
		@params:
			`gene`: The query gene
		@returns:
			A generator of `Enrichr_Library`
		"""
		response = requests.get(Enrichr.URL_GENEMAP.format(gene = gene))
		if not response.ok:
			raise Exception(Enrichr.MSG_ERROR_SEARCH_TERM)

		data = json.loads(response.text)
		descriptions = {desc['name']:(desc['description'] if 'description' in desc else '') \
			for desc in data['descriptions']}
		for category in data['categories']:
			for lib in category['libraries']:
				if lib['name'] == self.library:
					return Enrichr_Library(
						name        = lib['name'],
						category    = category['name'],
						hasGrid     = lib['hasGrid'],
						isFuzzy     = lib['isFuzzy'],
						format      = lib['format'],
						description = descriptions[self.library] \
							if self.library in descriptions else '',
						terms       = data['gene'][self.library] \
							if self.library in data['gene'] else []
					)

	def export(self, outfile, top = None):

		with open(outfile, 'w') as fout:
			fout.write('\t'.join(Enrichr_Term._fields) + '\n')
			if not self.results:
				return

			top = top or self.top
			results = self.results if top is None else self.results[:top]
			results = [ret for ret in results if getattr(ret, self.cutoff['by']) < self.cutoff['value']]

			for term in results:
				fout.write('\t'.join([
					"%.2E" % getattr(term, field) if 'Pval' in field else \
					"%.3f" % getattr(term, field) if field in ['Z', 'CombinedScore'] else \
					str(getattr(term, field)) for field in term._fields
				]) + '\n')

	def plot(self, outfile, top = None, res = 300, width = 2000, height = 2000):
		retfile = NamedTemporaryFile(delete = False)
		retfile.close()
		self.export(retfile.name, top)

		cmd = [
			self.Rscript,
			path.join(path.dirname(path.realpath(__file__)), 'plot.R'),
			retfile.name,
			outfile,
			str(self.size),
			str(res),
			str(width),
			str(height)
		]
		check_output(cmd)

		# import matplotlib, math
		# matplotlib.use('Agg')
		# from matplotlib import rcParams
		# rcParams['font.family'] = 'Arial'
		# from matplotlib import pyplot as plt
		# from matplotlib import gridspec
		# from matplotlib import patches

		# # already sorted
		# results = self.results[:top]

		# gs = gridspec.GridSpec(1, 2, width_ratios=[3, 7])
		# rownames = [r.term if len(r.term)<40 else r.term[:36] + ' ...' for r in results]
		# rnidx    = range(len(rownames))
		# ax1 = plt.subplot(gs[0])
		# plt.title (title.format(library = self.library), fontweight='bold')

		# ax1.xaxis.grid(alpha=.6, ls = '--', zorder = -99)
		# plt.subplots_adjust(wspace=.01, left=0.5)
		# ax1.barh(rnidx, [len(r.overlappingGenes) for r in results], color='blue', alpha=.6)
		# plt.yticks (rnidx, rownames)
		# ax1.yaxis.set_ticks_position('none')
		# ax1.tick_params(axis='x', colors='blue')
		# ax1.spines['top'].set_visible(False)
		# ax1.spines['left'].set_visible(False)
		# ax1.spines['right'].set_visible(False)
		# ax1.spines['bottom'].set_linewidth(1)
		# ax1.invert_xaxis()
		# ax1.invert_yaxis()
		# xticks = ax1.xaxis.get_major_ticks()
		# xticks[0].label1.set_visible(False)

		# ax2 = plt.subplot(gs[1])
		# ax2.xaxis.grid(alpha=.6, ls = '--', zorder = -99)
		# ax2.barh(rnidx, [-math.log(r[2], 10) for r in results], color='red', alpha = .6)
		# for i, r in enumerate(results):
		# 	t  = str("%.2E" % r[2])
		# 	tx = 0.1
		# 	ty = i + 0.1
		# 	ax2.text(tx, ty, t, fontsize=8)
		# ax2.tick_params(axis='x', colors='red')
		# ax2.spines['top'].set_visible(False)
		# ax2.spines['left'].set_visible(False)
		# ax2.spines['right'].set_visible(False)
		# ax2.spines['bottom'].set_linewidth(1)
		# ax2.invert_yaxis()
		# plt.yticks([])
		# ng_patch = patches.Patch(color='blue', alpha=.6, label='# overlapping genes')
		# pv_patch = patches.Patch(color='red', alpha = .6, label='-log(p-value)')
		# plt.figlegend(handles=[ng_patch, pv_patch], labels=['# overlapping genes', '-log(p-value)'], loc="lower center", ncol=2, edgecolor="none")
		# plt.savefig(outfile, dpi=300)
