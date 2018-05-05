import testly, tempfile, itertools
from os import path
from enrichr import Enrichr, Enrichr_Library, Enrichr_Term

GENES = [
	'PHF14', 'RBM3', 'MSL1', 'PHF21A', 'ARL10', 'INSR', 'JADE2', 'P2RX7',
	'LINC00662', 'CCDC101', 'PPM1B', 'KANSL1L', 'CRYZL1', 'ANAPC16', 'TMCC1',
	'CDH8', 'RBM11', 'CNPY2', 'HSPA1L', 'CUL2', 'PLBD2', 'LARP7', 'TECPR2',
	'ZNF302', 'CUX1', 'MOB2', 'CYTH2', 'SEC22C', 'EIF4E3', 'ROBO2',
	'ADAMTS9-AS2', 'CXXC1', 'LINC01314', 'ATF7', 'ATP5F1'
]
GENES_FILE   = tempfile.NamedTemporaryFile(delete = False)
RESULTS_FILE = tempfile.NamedTemporaryFile(delete = False)
RESULTS = [
	Enrichr_Term(rank=1, term='Longevity regulating pathway - multiple species_Homo sapiens_hsa04213', pval=0.00560, zscore=-2.01796, combinedScore=10.46188, overlappingGenes=['HSPA1L', 'INSR'], adjPval=0.22414, oldPval=0.0023517, oldAdjPval=0.09407),
	Enrichr_Term(rank=2, term='HIF-1 signaling pathway_Homo sapiens_hsa04066', pval=0.01399, zscore=-1.81648, combinedScore=7.75505, overlappingGenes=['INSR', 'CUL2'], adjPval=0.2798406, oldPval=0.005854, oldAdjPval=0.117082),
	Enrichr_Term(rank=3, term='Phospholipase D signaling pathway_Homo sapiens_hsa04072', pval=0.026229, zscore=-1.84575, combinedScore=6.720149, overlappingGenes=['CYTH2', 'INSR'], adjPval=0.33931, oldPval=0.0110856, oldAdjPval=0.1478086)
]

class TestEnrichr(testly.TestCase):

	def dataProvider_testInit(self):
		yield 'KEGG_2015',
		yield 'KEGG_2016',

	def testInit(self, library):
		enrichr = Enrichr(library)
		self.assertIsInstance(enrichr, Enrichr)
		self.assertIsNone(enrichr.listid)
		self.assertEqual(enrichr.library, library)
		self.assertEqual(enrichr.results, [])

	def dataProvider_testAddList(self):
		yield GENES, 'description'

	def testAddList(self, genes, description):
		enrichr = Enrichr()
		ret = enrichr.addList(genes, description)
		self.assertIsInstance(ret.userListId, int)
		self.assertEqual(ret.userListId, enrichr.listid)

	def dataProvider_testAddListFromFile(self):
		yield testly.Data(GENES_FILE.name, col = 1)

	def testAddListFromFile(self, genefile, delimit = '\t', skip = 0, col = 0):
		enrichr = Enrichr()
		ret = enrichr.addListFromFile(genefile, col = col, delimit = delimit, skip = skip)
		self.assertIsInstance(ret.userListId, int)
		self.assertEqual(ret.userListId, enrichr.listid)

	def dataProvider_testView(self):
		yield GENES, 'Another description'

	def testView(self, genes, description):
		self.maxDiff = None
		enrichr = Enrichr()
		enrichr.addList(genes, description)
		ret = enrichr.view()
		self.assertCountEqual(ret.genes, genes)
		self.assertEqual(ret.description, description)

	def dataProvider_testEnrich(self):
		yield GENES, 'Description', RESULTS

	def testEnrich(self, genes, description, results):
		self.maxDiff = None
		enrichr = Enrichr()
		enrichr.addList(genes, description)
		ret = enrichr.enrich()
		for i, result in enumerate(results):
			for field in result._fields:
				rval = getattr(result, field)
				oval = getattr(ret[i], field)
				if isinstance(rval, float):
					self.assertLess(abs(rval - oval), .01)
				else:
					self.assertEqual(rval, oval)

	def dataProvider_testGenemap(self):
		yield 'AKT1', [
			Enrichr_Library(name='ChEA_2016', category=u'Transcription', hasGrid=True, isFuzzy=True, format='{1} binds to the promoter region of {0}.', description='', terms=[
				'EGR1_19374776_ChIP-ChIP_THP-1_Human', 'CLOCK_20551151_ChIP-Seq_293T_Human',
				'SMARCA4_20176728_ChIP-ChIP_TSCs_Mouse', 'HNF4A_19761587_ChIP-ChIP_CACO-2_Human',
				'E2F1_21310950_ChIP-Seq_MCF-7_Human', 'EKLF_21900194_ChIP-Seq_ERYTHROCYTE_Mouse',
				'TCF3_18467660_ChIP-ChIP_MESCs_Mouse', 'ESRRB_18555785_ChIP-Seq_MESCs_Mouse',
				'E2F1_17053090_ChIP-ChIP_MCF-7_Human', 'GATA1_19941827_ChIP-Seq_MEL_Mouse',
				'TAL1_20566737_ChIP-Seq_PRIMARY_FETAL_LIVER_ERYTHROID_Mouse',
				'OLIG2_26023283_ChIP-Seq_AINV15_Mouse', 'FOXM1_26456572_ChIP-Seq_MCF-7_Human',
				'SMC4_20622854_ChIP-Seq_HELA_Human', 'DPY_21335234_ChIP-Seq_ESCs_Mouse',
				'TBX20_22080862_ChIP-Seq_HEART_Mouse', 'SOX3_22085726_ChIP-Seq_MUSCLE_Mouse',
				'RXR_22158963_ChIP-Seq_LIVER_Mouse', 'TBX20_22328084_ChIP-Seq_HEART_Mouse',
				'GATA1_22383799_ChIP-Seq_G1ME_Mouse', 'SMC1_22415368_ChIP-Seq_MEFs_Mouse',
				'SMC3_22415368_ChIP-Seq_MEFs_Mouse', 'NCOR_22465074_ChIP-Seq_MACROPHAGES_Mouse',
				'CTCF_26484167_Chip-Seq_Bcells_Mouse', 'MAF_26560356_Chip-Seq_TH2_Human',
				'CREB1_26743006_Chip-Seq_LNCaP-abl_Human', 'CREB1_26743006_Chip-Seq_LNCaP_Human',
				'KLF4_26769127_Chip-Seq_PDAC-Cell_line_Human', 'KDM2B_26808549_Chip-Seq_HPB-ALL_Human',
				'KDM2B_26808549_Chip-Seq_SIL-ALL_Human', 'RACK7_27058665_Chip-Seq_MCF-7_Human',
				'BRD4_27068464_Chip-Seq_AML-cells_Mouse', 'ATF3_27146783_Chip-Seq_COLON_Human',
				'CTCF_27219007_Chip-Seq_Bcells_Human', 'CTCF_27219007_Chip-Seq_ERYTHROID_Human',
				'KAP1_27257070_Chip-Seq_ESCs_Mouse', 'SOX2_27498859_Chip-Seq_STOMACH_Mouse',
				'EP300_20729851_ChIP-Seq_FORBRAIN_MIDBRAIN_LIMB_HEART_Mouse',
				'POU5F1_18347094_ChIP-ChIP_MESCs_Mouse', 'REST_19997604_ChIP-ChIP_NEURONS_Mouse',
				'TBX5_21415370_ChIP-Seq_HL-1_Mouse', 'YAP1_20516196_ChIP-Seq_MESCs_Mouse',
				'RCOR1_19997604_ChIP-ChIP_NEURONS_Mouse',
				'TCFAP2C_20176728_ChIP-ChIP_TROPHOBLAST_STEM_CELLS_Mouse',
				'MTF2_20144788_ChIP-Seq_MESCs_Mouse',
				'CUX1_19635798_ChIP-ChIP_MULTIPLE_HUMAN_CANCER_TYPES_Human',
				'TP53_23651856_ChIP-Seq_MEFs_Mouse', 'ZFX_18555785_ChIP-Seq_MESCs_Mouse',
				'OCT4_18692474_ChIP-Seq_MEFs_Mouse', 'POU5F1_18692474_ChIP-Seq_MESCs_Mouse',
				'TP63_22573176_ChIP-Seq_HFKS_Human', 'SUZ12_20075857_ChIP-Seq_MESCs_Mouse',
				'HNF4A_19822575_ChIP-Seq_HepG2_Human', 'EGR1_20690147_ChIP-Seq_ERYTHROLEUKEMIA_Human']),
			Enrichr_Library(name='TRANSFAC_and_JASPAR_PWMs', category='Transcription', hasGrid=True, isFuzzy=False, format='{1} has a binding site at the promoter of {0}.', description='PWMs from TRANSFAC and JASPAR were used to scan the promoters of all human genes in the region ?2000 and +500 from the transcription factor start site (TSS). We retained only the 100% matches to the consensus sequences to call an interaction between a factor and target gene. This gene-set library was created for a tool we previously published called Expression2Kinases (PMID: 22080467).', terms=[
				'EGR1 (mouse)', 'TEAD1 (human)', 'CACYBP (human)', 'SMAD4 (mouse)', 'ELF3 (human)', 'SP3 (human)', 'TFAP2A (human)', 'SP1 (mouse)', 'ETV4 (human)', 'NR5A2 (mouse)',
				'ZNF148 (human)', 'RUNX1 (human)', 'MIR138 (human)', 'GATA2 (human)', 'STAT5B (mouse)', 'RELA (human)', 'TCFAP2A (human)', 'NFE2 (human)', 'RBPJ (human)',
				'PCBP1 (human)', 'HINFP (human)', 'YY1 (mouse)', 'MZF1_1-4 (human)', 'ETS1 (human)', 'AHR (human)', 'CRTC1 (human)', 'SP1 (human)', 'NFKB1 (human)']),
			Enrichr_Library(name='ARCHS4_TFs_Coexp', category='Transcription', hasGrid=False, isFuzzy=True, format='{0} is coexpressed with {1}.', description='', terms=[
				'TEAD3_human_tf_ARCHS4_coexpression', 'HOXC8_human_tf_ARCHS4_coexpression',
				'YBX1_human_tf_ARCHS4_coexpression', 'TFE3_human_tf_ARCHS4_coexpression',
				'TCF3_human_tf_ARCHS4_coexpression', 'SREBF1_human_tf_ARCHS4_coexpression',
				'CIC_human_tf_ARCHS4_coexpression', 'NFE2L1_human_tf_ARCHS4_coexpression',
				'GRHL2_human_tf_ARCHS4_coexpression', 'RCOR1_human_tf_ARCHS4_coexpression',
				'ESR1_human_tf_ARCHS4_coexpression', 'CEBPE_human_tf_ARCHS4_coexpression',
				'MZF1_human_tf_ARCHS4_coexpression', 'XBP1_human_tf_ARCHS4_coexpression',
				'TCFL5_human_tf_ARCHS4_coexpression', 'GMEB2_human_tf_ARCHS4_coexpression',
				'E2F1_human_tf_ARCHS4_coexpression', 'MAZ_human_tf_ARCHS4_coexpression',
				'USF2_human_tf_ARCHS4_coexpression', 'ERF_human_tf_ARCHS4_coexpression',
				'GATA3_human_tf_ARCHS4_coexpression', 'UBTF_human_tf_ARCHS4_coexpression',
				'WHSC1_human_tf_ARCHS4_coexpression', 'PRDM4_human_tf_ARCHS4_coexpression',
				'SRF_human_tf_ARCHS4_coexpression', 'HOXB3_human_tf_ARCHS4_coexpression',
				'MSX2_human_tf_ARCHS4_coexpression', 'ZNF205_human_tf_ARCHS4_coexpression',
				'HOXC13_human_tf_ARCHS4_coexpression', 'HOXC11_human_tf_ARCHS4_coexpression',
				'MXD4_human_tf_ARCHS4_coexpression', 'TRERF1_human_tf_ARCHS4_coexpression',
				'SPDEF_human_tf_ARCHS4_coexpression', 'ZNF384_human_tf_ARCHS4_coexpression',
				'ELK1_human_tf_ARCHS4_coexpression', 'ATF4_human_tf_ARCHS4_coexpression',
				'FOXA1_human_tf_ARCHS4_coexpression', 'KLF16_human_tf_ARCHS4_coexpression',
				'TFAP2A_human_tf_ARCHS4_coexpression', 'CREBZF_human_tf_ARCHS4_coexpression',
				'ESR2_human_tf_ARCHS4_coexpression', 'FOXK2_human_tf_ARCHS4_coexpression',
				'RFX5_human_tf_ARCHS4_coexpression', 'UBP1_human_tf_ARCHS4_coexpression',
				'SKI_human_tf_ARCHS4_coexpression', 'ZBTB7B_human_tf_ARCHS4_coexpression',
				'SNAPC4_human_tf_ARCHS4_coexpression', 'SMAD3_human_tf_ARCHS4_coexpression',
				'E4F1_human_tf_ARCHS4_coexpression', 'FOXI1_human_tf_ARCHS4_coexpression',
				'ZNF282_human_tf_ARCHS4_coexpression', 'HOXC5_human_tf_ARCHS4_coexpression',
				'RARG_human_tf_ARCHS4_coexpression', 'RELA_human_tf_ARCHS4_coexpression',
				'ESRRA_human_tf_ARCHS4_coexpression', 'SIX5_human_tf_ARCHS4_coexpression',
				'IRX3_human_tf_ARCHS4_coexpression', 'HSF1_human_tf_ARCHS4_coexpression',
				'SP1_human_tf_ARCHS4_coexpression', 'RXRA_human_tf_ARCHS4_coexpression',
				'ZNF395_human_tf_ARCHS4_coexpression', 'HOXC6_human_tf_ARCHS4_coexpression',
				'HOXC4_human_tf_ARCHS4_coexpression', 'MAFK_human_tf_ARCHS4_coexpression',
				'SREBF2_human_tf_ARCHS4_coexpression', 'RXRB_human_tf_ARCHS4_coexpression',
				'PBX2_human_tf_ARCHS4_coexpression', 'E2F4_human_tf_ARCHS4_coexpression',
				'ZBED1_human_tf_ARCHS4_coexpression', 'MEOX1_human_tf_ARCHS4_coexpression',
				'ZNF263_human_tf_ARCHS4_coexpression', 'WIZ_human_tf_ARCHS4_coexpression',
				'ZNF629_human_tf_ARCHS4_coexpression', 'GLIS2_human_tf_ARCHS4_coexpression',
				'MLLT1_human_tf_ARCHS4_coexpression', 'FOXP4_human_tf_ARCHS4_coexpression',
				'ZNF740_human_tf_ARCHS4_coexpression', 'ZNF710_human_tf_ARCHS4_coexpression',
				'ZNF687_human_tf_ARCHS4_coexpression', 'CIZ1_human_tf_ARCHS4_coexpression',
				'PRDM8_human_tf_ARCHS4_coexpression', 'CUX1_human_tf_ARCHS4_coexpression',
				'ZNF653_human_tf_ARCHS4_coexpression', 'ZNF592_human_tf_ARCHS4_coexpression',
				'ZNF217_human_tf_ARCHS4_coexpression', 'IRX5_human_tf_ARCHS4_coexpression',
				'ZNF664_human_tf_ARCHS4_coexpression', 'FIZ1_human_tf_ARCHS4_coexpression',
				'ZNF609_human_tf_ARCHS4_coexpression', 'ZFP41_human_tf_ARCHS4_coexpression',
				'ZNF703_human_tf_ARCHS4_coexpression', 'ZFP91_human_tf_ARCHS4_coexpression',
				'ZNF777_human_tf_ARCHS4_coexpression', 'ZBTB22_human_tf_ARCHS4_coexpression',
				'ZBTB12_human_tf_ARCHS4_coexpression', 'TCF20_human_tf_ARCHS4_coexpression',
				'TCF25_human_tf_ARCHS4_coexpression', 'ZC3H3_human_tf_ARCHS4_coexpression',
				'NCOA3_human_tf_ARCHS4_coexpression', 'ZNF787_human_tf_ARCHS4_coexpression',
				'ZMAT2_human_tf_ARCHS4_coexpression', 'FOXK1_human_tf_ARCHS4_coexpression',
				'ZNF598_human_tf_ARCHS4_coexpression', 'ZNF768_human_tf_ARCHS4_coexpression',
				'ZNF316_human_tf_ARCHS4_coexpression', 'PLXND1_human_tf_ARCHS4_coexpression',
				'DVL2_human_tf_ARCHS4_coexpression', 'CSDE1_human_tf_ARCHS4_coexpression',
				'AKAP8L_human_tf_ARCHS4_coexpression', 'BRD9_human_tf_ARCHS4_coexpression',
				'HMG20B_human_tf_ARCHS4_coexpression', 'MBD3_human_tf_ARCHS4_coexpression',
				'SRCAP_human_tf_ARCHS4_coexpression', 'ZC3H7B_human_tf_ARCHS4_coexpression',
				'DOT1L_human_tf_ARCHS4_coexpression', 'SF3A2_human_tf_ARCHS4_coexpression',
				'DVL1_human_tf_ARCHS4_coexpression', 'PLXNA1_human_tf_ARCHS4_coexpression',
				'SLC2A4RG_human_tf_ARCHS4_coexpression', 'PRR12_human_tf_ARCHS4_coexpression',
				'UBR4_human_tf_ARCHS4_coexpression', 'PLXNA3_human_tf_ARCHS4_coexpression',
				'UNK_human_tf_ARCHS4_coexpression', 'SMAD6_human_tf_ARCHS4_coexpression',
				'SMARCC2_human_tf_ARCHS4_coexpression', 'MINK1_human_tf_ARCHS4_coexpression',
				'MBD1_human_tf_ARCHS4_coexpression', 'POGK_human_tf_ARCHS4_coexpression',
				'SETDB1_human_tf_ARCHS4_coexpression', 'MTA2_human_tf_ARCHS4_coexpression',
				'BRPF1_human_tf_ARCHS4_coexpression', 'ADAR_human_tf_ARCHS4_coexpression',
				'DVL3_human_tf_ARCHS4_coexpression', 'PLXNB1_human_tf_ARCHS4_coexpression',
				'ATMIN_human_tf_ARCHS4_coexpression', 'GATAD2A_human_tf_ARCHS4_coexpression',
				'MECP2_human_tf_ARCHS4_coexpression', 'SSH3_human_tf_ARCHS4_coexpression',
				'KAT5_human_tf_ARCHS4_coexpression', 'ANAPC2_human_tf_ARCHS4_coexpression',
				'CCDC71_human_tf_ARCHS4_coexpression', 'TIGD5_human_tf_ARCHS4_coexpression',
				'RBM10_human_tf_ARCHS4_coexpression', 'MTA1_human_tf_ARCHS4_coexpression',
				'EP400_human_tf_ARCHS4_coexpression', 'H1FX_human_tf_ARCHS4_coexpression',
				'H1F0_human_tf_ARCHS4_coexpression', 'NCOR2_human_tf_ARCHS4_coexpression',
				'PLXNB2_human_tf_ARCHS4_coexpression', 'ZNF512B_human_tf_ARCHS4_coexpression',
				'INF2_human_tf_ARCHS4_coexpression', 'REPIN1_human_tf_ARCHS4_coexpression'])]

	def testGenemap(self, gene, results):
		enrichr = Enrichr()
		ret = list(itertools.islice(enrichr.genemap(gene), 3))
		self.assertListEqual(ret, results)

	def dataProvider_testExport(self):
		yield GENES, 'Description', RESULTS_FILE.name, RESULTS

	def testExport(self, genes, description, retfile, results):
		self.maxDiff = None
		enrichr = Enrichr()
		enrichr.addList(genes, description)
		enrichr.enrich()
		enrichr.export(retfile)
		ret = []
		with open(retfile) as f:
			for line in f:
				line = line.strip()
				if not line or line.startswith('rank'): continue
				ret.append(Enrichr_Term._make(line.split('\t')))
		for i, result in enumerate(results):
			for field in result._fields:
				rval = getattr(result, field)
				oval = getattr(ret[i], field)
				if field == 'overlappingGenes':
					oval = oval.split('|')
				if isinstance(rval, (int, float)):
					self.assertLess(abs(rval - float(oval)), .01)
				else:
					self.assertEqual(rval, oval)

	def dataProvider_testPlot(self):
		yield GENES, 'Description', path.join(tempfile.gettempdir(), 'plot.png')

	def testPlot(self, genes, description, plotfile):
		enrichr = Enrichr()
		enrichr.addList(genes, description)
		enrichr.enrich()
		enrichr.plot(plotfile)
		self.assertGreater(path.getsize(plotfile), 0)

if __name__ == '__main__':
	GENES_FILE.write('\n'.join(['%s\t%s' % (i, gene) for i, gene in enumerate(GENES)]).encode())
	GENES_FILE.close()
	testly.main(verbosity = 2)
