"""API for Enrichr"""

import re
import sys
import json
from tempfile import NamedTemporaryFile
from collections import namedtuple
from subprocess import check_output
from os import path
import requests
VERSION = '0.0.0dev'


# types
EnrichrListId = namedtuple('EnrichrListId', ['userListId', 'shortId'])
EnrichrGeneList = namedtuple('EnrichrGeneList', ['genes', 'description'])
EnrichrTerm = namedtuple('EnrichrTerm', [
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
EnrichrLibrary = namedtuple('EnrichrLibrary', [
    'name',
    'category',
    'hasGrid',
    'isFuzzy',
    'format',
    'description',
    'terms'
])

# urls
URL_ADDLIST = 'http://amp.pharm.mssm.edu/Enrichr/addList'
URL_VIEW = 'http://amp.pharm.mssm.edu/Enrichr/view?userListId={listid}'
URL_ENRICH = ('http://amp.pharm.mssm.edu/Enrichr/enrich?userListId={listid}'
              '&backgroundType={library}')
URL_GENEMAP = ('http://amp.pharm.mssm.edu/Enrichr/genemap?json=true&'
               'setup=true&gene={gene}')
URL_EXPORT = ('http://amp.pharm.mssm.edu/Enrichr/export?userListId={listid}'
              '&backgroundType={library}')

# msgs
MSG_NO_LISTID = 'You have to successfully add a list first'
MSG_NO_LIBRARY = 'You have to set a library'
MSG_NO_RESULTS = 'You have to run an enrich first'
MSG_ERROR_ADD_GENELIST = 'Error adding gene list'
MSG_ERROR_GET_GENELIST = 'Error getting gene list'
MSG_ERROR_ENRICH = 'Error fetching enrichment results'
MSG_ERROR_SEARCH_TERM = 'Error searching for terms'

class Enrichr:

    """Main class"""

    def __init__(self, library='KEGG_2016',
                 cutoff=0.05, top=20, rscript='Rscript'):
        """
        Constructor
        @params:
            `library`: The library to do enrichment against.
                Default: `KEGG_2016`
                - See all available libraries:
                  http://amp.pharm.mssm.edu/Enrichr/#stats
            `cutoff` : The cutoff for the enriched terms. Default: `0.05`
                - A float indicating cutoff by "AdjPval", or
                - A dict specifying the `value` and `by`
            `top`: Top N terms. Filtering terms with `cutoff`,
                   whichever comes first. Default: `20`
                - Filters (cutoff and top) are for export and plot
        """
        # save the session
        self.listid = None
        self.size = None
        self.library = library
        self.results = []
        self.cutoff = cutoff if isinstance(cutoff, dict) else {
            "by": "AdjPval", "value": cutoff}
        self.top = top
        self.rscript = rscript

    def addListFromFile(self, genefile, # pylint: disable=invalid-name
                        col=0, delimit='\t', skip=0):
        """
        Add a gene list from files
        @params:
            `genefile`: The file containing genes
            `col`     : The column index where the genes located. Default: `0`
            `delimit` : The delimiter for the columns. Default : `\t`
            `skip`    : Skip first N lines. Default : `0`
        @returns:
            A `EnrichrListId` with the session ids.
        """
        description = path.basename(genefile)
        genes = set()
        with open(genefile) as fgene:
            for _ in range(skip):
                fgene.readline()
            for line in fgene:
                line = line.strip()
                if not line:
                    continue
                genes.add(line.split(delimit)[col].strip())
        return self.addList(genes, description)
    add_list_from_file = addListFromFile

    def addList(self, genes,  # pylint: disable=invalid-name
                description='Enrichr gene list.'):
        """
        Add a gene list
        @params:
            `genes`: The genes.
            `description`: The description of the gene list
        @returns:
            A `EnrichrListId` with the session ids.
        """
        self.size = len(genes)
        genes_str = '\n'.join(set(genes))
        payload = {
            'list': (None, genes_str),
            'description': (None, description)
        }
        response = requests.post(URL_ADDLIST, files=payload)
        if not response.ok:
            msg = re.sub('<[^<]+?>', '', response.text).splitlines()
            #msg = response.text.splitlines()
            msg = [line.strip() for line in msg if line.strip()]
            msg = ', '.join(msg)
            raise Exception(MSG_ERROR_ADD_GENELIST + ': ' + msg)

        data = json.loads(response.text)
        self.listid = data['userListId']
        return EnrichrListId(**data)
    add_list = addList

    def view(self, listid=None):
        """
        View added gene list
        @params:
            `listid`: The id of the gene list. Default: `None` (current list)
        @returns:
            `EnrichrGeneList`
        """
        listid = listid or self.listid
        if not listid:
            raise Exception(MSG_NO_LISTID)
        self.listid = listid

        response = requests.get(URL_VIEW.format(listid=listid))
        if not response.ok:
            raise Exception(MSG_ERROR_GET_GENELIST)

        data = json.loads(response.text)
        return EnrichrGeneList(genes=data['genes'],
                               description=data['description'])

    def enrich(self, library=None, listid=None):
        """
        Do the enrichment
        @params:
            `library`: THe library to do enrichment against.
                       Default: `None` (current library)
            `listid` : The list id of the session.
                       Default: `None` (current session)
        @returns:
            A list of `EnrichrTerm` (sorted by `adjPval`)
        """
        listid = listid or self.listid
        if not listid:
            raise Exception(MSG_NO_LISTID)
        self.listid = listid

        library = library or self.library
        if not library:
            raise Exception(MSG_NO_LIBRARY)
        self.library = library

        response = requests.get(URL_EXPORT.format(
            listid=listid, library=library), stream=True)
        response.encoding = 'utf-8'
        if "Error" in response.text:
            raise RuntimeError('%s\n%s' % (
                self.library, re.sub(r'\s*<[^<]+?>\s*', ' ', response.text)
            ))

        self.results = []
        for line in response.text.splitlines():
            if line.startswith('Term'):
                continue
            data = line.split('\t')
            data[2:-1] = [float(d) for d in data[2:-1]]
            data[0] = data[0].replace('_Homo sapiens', '')
            self.results.append(EnrichrTerm._make(data))
        self.results = sorted(self.results, key=lambda x: x.AdjPval)
        return self.results

    def genemap(self, gene, libs=None):
        """
        Find terms that contain a given gene
        @params:
            `gene`: The query gene
        @returns:
            A generator of `EnrichrLibrary`
        """
        response = requests.get(URL_GENEMAP.format(gene=gene))
        if not response.ok:
            raise Exception(MSG_ERROR_SEARCH_TERM)
        lib_query = libs or self.library
        if not isinstance(lib_query, list):
            lib_query = [lib_query]

        data = json.loads(response.text)
        descriptions = {desc['name']: (desc['description']
                                       if 'description' in desc else '')
                        for desc in data['descriptions']}
        ret = {}
        for category in data['categories']:
            for lib in category['libraries']:
                libname = lib['name']
                if not libname in lib_query:
                    continue

                ret[libname] = EnrichrLibrary(
                    name=libname,
                    category=category['name'],
                    hasGrid=lib['hasGrid'],
                    isFuzzy=lib['isFuzzy'],
                    format=lib['format'],
                    description=descriptions[libname]
                    if libname in descriptions else '',
                    terms=data['gene'][libname]
                    if libname in data['gene'] else []
                )

        return (ret if libs else ret[self.library]
                if self.library in ret else None)

    def export(self, outfile, top=None):
        """Export the enriched pathways"""

        with open(outfile, 'w') as fout:
            fout.write('\t'.join(EnrichrTerm._fields) + '\n')
            if not self.results:
                return

            top = top or self.top
            results = self.results if top is None else self.results[:top]
            results = [ret for ret in results if getattr(
                ret, self.cutoff['by']) < self.cutoff['value']]

            for term in results:
                fout.write('\t'.join([
                    "%.2E" % getattr(term, field)
                    if 'Pval' in field
                    else "%.3f" % getattr(term, field)
                    if field in ['Z', 'CombinedScore']
                    else str(getattr(term, field))
                    for field in term._fields
                ]) + '\n')

    def plot(self,  # pylint: disable=too-many-arguments
             outfile,
             top=None,
             res=300,
             width=2000,
             height=2000):
        """Plot the enrichment figure"""
        retfile = NamedTemporaryFile(delete=False)
        retfile.close()
        self.export(retfile.name, top)

        cmd = [
            self.rscript,
            path.join(path.dirname(path.realpath(__file__)), 'plot.R'),
            retfile.name,
            outfile,
            str(self.size),
            str(res),
            str(width),
            str(height)
        ]
        sys.stderr.write("[INFO] Generating plots by running: \n")
        sys.stderr.write("[INFO]   %s\n" % ' '.join(cmd))
        check_output(cmd)
