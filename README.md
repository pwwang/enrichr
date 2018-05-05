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
>>> results = en.enrich()
>>> # perform enrichment on a different library and a different gene list:
>>> # results = en.enrich(library = 'ChEA_2016', listid = 123456)
>>> print results # a list of namedtuple
[Enrichr_Term(rank=1, term=u'Longevity regulating pathway - multiple species_Homo sapiens_hsa04213', pval=0.00560341272708137, zscore=-2.017962732978835, combinedScore=10.461884526362898, overlappingGenes=[u'HSPA1L', u'INSR'], adjPval=0.2241365090832548, oldPval=0.0023516996733944644, oldAdjPval=0.09406798693577857), Enrichr_Term(...), ...]
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
rank	term	pval	zscore	combinedScore	overlappingGenes	adjPval	oldPval	oldAdjPval
1	Longevity regulating pathway - multiple species_Homo sapiens_hsa04213	0.00560341272708	-2.01796273298	10.4618845264	HSPA1L|INSR	0.224136509083	0.00235169967339	0.0940679869358
2	HIF-1 signaling pathway_Homo sapiens_hsa04066	0.0139920292251	-1.81648257202	7.75504992289	INSR|CUL2	0.279840584501	0.00585409985605	0.117081997121
3	Phospholipase D signaling pathway_Homo sapiens_hsa04072	0.0262293541331	-1.84575051838	6.72014896483	CYTH2|INSR    0.339310255604	0.0110856438446	0.147808584594
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
