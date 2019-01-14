# enrichr: A python wrapper for Enrichr APIs
[Enrichr][1] is a web tool for gene set enrichment analysis.

## Features
- Compatible with both python2 and python3
- Able to export part of the enrichment results
- Able to plot the enrichment results

## Installation
```shell
> pip install git+https://github.com/pwwang/enrichr.git
> # will install requests, matplotlib
> # to run tests, you have to install testly
> # pip install git+https://github.com/pwwang/testly.git
> python test.py
```

## Usage
### Initiate an `enrichr` instance
```python
>>> from enrichr import Enrichr
>>> en = Enrichr(library = 'KEGG_2016')
>>> # use a different library:
>>> # en = Enrichr('ChEA_2016')
>>> # Add different cutoff
>>> en = Enrichr(cutoff = 0.01, top = 10)
>>> # Use adjusted p-value 0.01 
>>> # whichever comes first, records with adjusted p-value or top
```

### Add a gene list
```python
>>> en.addList(['PHF14', 'RBM3', 'MSL1', ...])
>>> # set a description and show the userListId:
>>> # r = en.addList([...], description = 'Some description')
>>> # print r.userListId
>>> # 7260214
```

### Add a gene list from a file
```shell
> head genes.txt
#ID,Gene
1,PHF14
2,RBM3
3,MSL1
```
```python
>>> en.addListFromFile('genes.txt', col = 1, delimit = ',', skip = 1)
```

### View added gene list
```python
>>> v = en.view()
>>> # if you have a different listid:
>>> # v = en.view(listid = 1234567)
>>> print v # namedtuple
Enrichr_Genelist(genes=[u'PHF14', u'RBM3', u'MSL1', ...], description=u'Enrichr gene list.')
>>> print v.genes, v.description
[u'PHF14', u'RBM3', u'MSL1', ...], 'Enrichr gene list.'
```

### Do the enrichment
```python
>>> results = en.enrich(cutoff = .25)
>>> # perform enrichment on a different library and a different gene list:
>>> # results = en.enrich(library = 'ChEA_2016', listid = 123456)
>>> print results # a list of namedtuple
[Enrichr_Term(
	Term          = 'Longevity regulating pathway - multiple species_hsa04213',
	Overlap       = '2/64',
	Pval          = 0.00560341272708137,
	AdjPval       = 0.2241365090832548,
	OldPval       = 0.0023516996733944644,
	OldAdjPval    = 0.09406798693577857,
	Z             = -2.017962732978835,
	CombinedScore = 10.461884526362898,
	Genes         = 'HSPA1L;INSR'), Enrichr_Term(...), ...]
```

### Export the results
```python
>>> en.enrich()
>>> en.export('results.txt')
>>> # only export top 3 terms:
>>> # en.export('results.txt', top = 3)
```
```shell
> cat results.txt
Term	Overlap	Pval	AdjPval	OldPval	OldAdjPval	Z	CombinedScore	Genes
Longevity regulating pathway - multiple species_hsa04213	2/64	5.60E-03	2.24E-01	2.35E-03	9.41E-02	-2.018	10.462	HSPA1L;INSR
HIF-1 signaling pathway_hsa04066	2/103	1.40E-02	2.80E-01	5.85E-03	1.17E-01	-1.816	7.755	INSR;CUL2
Phospholipase D signaling pathway_hsa04072	2/144	2.62E-02	3.39E-01	1.11E-02	1.48E-01	-1.846	6.720	CYTH2;INSR
MAPK signaling pathway_hsa04010	2/255	7.32E-02	3.39E-01	3.23E-02	2.34E-01	-1.896	4.957	PPM1B;HSPA1L
...
```

### Plot the results
```python
>>> en.enrich()
>>> en.plot('results.png')
>>> # set a different title and different number of terms to plot:
>>> # en.plot('results.png', title = 'Gene enrichment: {library}', top = 20)
```
![results.png][2]

### Find terms that contain a given gene
```python
>>> libs = en.genemap('AKT1') # it's a generator
>>> print libs
<generator object genemap at 0x2b127b472730>
>>> lib = next(libs)
>>> print lib # namedtuple
Enrichr_Library(name='ChEA_2016', category='Transcription', hasGrid=True, isFuzzy=True, format='{1} binds to the promoter region of {0}.', description='', terms=['EGR1_19374776_ChIP-ChIP_THP-1_Human', 'CLOCK_20551151_ChIP-Seq_293T_Human', ...])
```

[1]: http://amp.pharm.mssm.edu/Enrichr/
[2]: ./results.png
