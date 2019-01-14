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
	Enrichr_Term(
		Term          = 'Longevity regulating pathway - multiple species_hsa04213',
		Overlap       = '2/64',
		Pval          = 0.00560341272708137,
		AdjPval       = 0.2241365090832548,
		OldPval       = 0.0023516996733944644,
		OldAdjPval    = 0.09406798693577857,
		Z             = -2.017962732978835,
		CombinedScore = 10.461884526362898,
		Genes         = 'HSPA1L;INSR'),
	Enrichr_Term(
		Term          = 'HIF-1 signaling pathway_hsa04066',
		Overlap       = '2/103',
		Pval          = 0.013992029225050977,
		AdjPval       = 0.27984058450101956,
		OldPval       = 0.005854099856053235,
		OldAdjPval    = 0.11708199712106471,
		Z             = -1.8164825720197595,
		CombinedScore = 7.755049922886117,
		Genes         = 'INSR;CUL2'),
	Enrichr_Term(
		Term          = 'Phospholipase D signaling pathway_hsa04072',
		Overlap       = '2/144',
		Pval          = 0.026229354133096315,
		AdjPval       = 0.33931025560360656,
		OldPval       = 0.011085643844551444,
		OldAdjPval    = 0.14780858459401924,
		Z             = -1.845750518377502,
		CombinedScore = 6.720148964834603,
		Genes         = 'CYTH2;INSR'),
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
				if isinstance(rval, (int, float)):
					self.assertAlmostEqual(rval, oval)
				else:
					self.assertEqual(rval, oval)


	def dataProvider_testGenemap(self):
		yield 'AKT1', [
			Enrichr_Library(name='ChEA_2016', category=u'Transcription', hasGrid=True, isFuzzy=False, format='{1} binds to the promoter region of {0}.', description='', terms=[
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
			Enrichr_Library(name=u'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X', category=u'Transcription', hasGrid=False, isFuzzy=False, format=u'{1} binds to the promoter region of {0}.', description='', terms=[u'EGR1_CHEA', u'E2F1_CHEA']),
			Enrichr_Library(name=u'ARCHS4_TFs_Coexp', category=u'Transcription', hasGrid=False, isFuzzy=False, format=u'{0} is coexpressed with {1}.', description='', terms=[u'TEAD3_human_tf_ARCHS4_coexpression', u'HOXC8_human_tf_ARCHS4_coexpression', u'YBX1_human_tf_ARCHS4_coexpression', u'TFE3_human_tf_ARCHS4_coexpression', u'TCF3_human_tf_ARCHS4_coexpression', u'SREBF1_human_tf_ARCHS4_coexpression', u'CIC_human_tf_ARCHS4_coexpression', u'NFE2L1_human_tf_ARCHS4_coexpression', u'GRHL2_human_tf_ARCHS4_coexpression', u'RCOR1_human_tf_ARCHS4_coexpression', u'ESR1_human_tf_ARCHS4_coexpression', u'CEBPE_human_tf_ARCHS4_coexpression', u'MZF1_human_tf_ARCHS4_coexpression', u'XBP1_human_tf_ARCHS4_coexpression', u'TCFL5_human_tf_ARCHS4_coexpression', u'GMEB2_human_tf_ARCHS4_coexpression', u'E2F1_human_tf_ARCHS4_coexpression', u'MAZ_human_tf_ARCHS4_coexpression', u'USF2_human_tf_ARCHS4_coexpression', u'ERF_human_tf_ARCHS4_coexpression', u'GATA3_human_tf_ARCHS4_coexpression', u'UBTF_human_tf_ARCHS4_coexpression', u'WHSC1_human_tf_ARCHS4_coexpression', u'PRDM4_human_tf_ARCHS4_coexpression', u'SRF_human_tf_ARCHS4_coexpression', u'HOXB3_human_tf_ARCHS4_coexpression', u'MSX2_human_tf_ARCHS4_coexpression', u'ZNF205_human_tf_ARCHS4_coexpression', u'HOXC13_human_tf_ARCHS4_coexpression', u'HOXC11_human_tf_ARCHS4_coexpression', u'MXD4_human_tf_ARCHS4_coexpression', u'TRERF1_human_tf_ARCHS4_coexpression', u'SPDEF_human_tf_ARCHS4_coexpression', u'ZNF384_human_tf_ARCHS4_coexpression', u'ELK1_human_tf_ARCHS4_coexpression', u'ATF4_human_tf_ARCHS4_coexpression', u'FOXA1_human_tf_ARCHS4_coexpression', u'KLF16_human_tf_ARCHS4_coexpression', u'TFAP2A_human_tf_ARCHS4_coexpression', u'CREBZF_human_tf_ARCHS4_coexpression', u'ESR2_human_tf_ARCHS4_coexpression', u'FOXK2_human_tf_ARCHS4_coexpression', u'RFX5_human_tf_ARCHS4_coexpression', u'UBP1_human_tf_ARCHS4_coexpression', u'SKI_human_tf_ARCHS4_coexpression', u'ZBTB7B_human_tf_ARCHS4_coexpression', u'SNAPC4_human_tf_ARCHS4_coexpression', u'SMAD3_human_tf_ARCHS4_coexpression', u'E4F1_human_tf_ARCHS4_coexpression', u'FOXI1_human_tf_ARCHS4_coexpression', u'ZNF282_human_tf_ARCHS4_coexpression', u'HOXC5_human_tf_ARCHS4_coexpression', u'RARG_human_tf_ARCHS4_coexpression', u'RELA_human_tf_ARCHS4_coexpression', u'ESRRA_human_tf_ARCHS4_coexpression', u'SIX5_human_tf_ARCHS4_coexpression', u'IRX3_human_tf_ARCHS4_coexpression', u'HSF1_human_tf_ARCHS4_coexpression', u'SP1_human_tf_ARCHS4_coexpression', u'RXRA_human_tf_ARCHS4_coexpression', u'ZNF395_human_tf_ARCHS4_coexpression', u'HOXC6_human_tf_ARCHS4_coexpression', u'HOXC4_human_tf_ARCHS4_coexpression', u'MAFK_human_tf_ARCHS4_coexpression', u'SREBF2_human_tf_ARCHS4_coexpression', u'RXRB_human_tf_ARCHS4_coexpression', u'PBX2_human_tf_ARCHS4_coexpression', u'E2F4_human_tf_ARCHS4_coexpression', u'ZBED1_human_tf_ARCHS4_coexpression', u'MEOX1_human_tf_ARCHS4_coexpression', u'ZNF263_human_tf_ARCHS4_coexpression', u'WIZ_human_tf_ARCHS4_coexpression', u'ZNF629_human_tf_ARCHS4_coexpression', u'GLIS2_human_tf_ARCHS4_coexpression', u'MLLT1_human_tf_ARCHS4_coexpression', u'FOXP4_human_tf_ARCHS4_coexpression', u'ZNF740_human_tf_ARCHS4_coexpression', u'ZNF710_human_tf_ARCHS4_coexpression', u'ZNF687_human_tf_ARCHS4_coexpression', u'CIZ1_human_tf_ARCHS4_coexpression', u'PRDM8_human_tf_ARCHS4_coexpression', u'CUX1_human_tf_ARCHS4_coexpression', u'ZNF653_human_tf_ARCHS4_coexpression', u'ZNF592_human_tf_ARCHS4_coexpression', u'ZNF217_human_tf_ARCHS4_coexpression', u'IRX5_human_tf_ARCHS4_coexpression', u'ZNF664_human_tf_ARCHS4_coexpression', u'FIZ1_human_tf_ARCHS4_coexpression', u'ZNF609_human_tf_ARCHS4_coexpression', u'ZFP41_human_tf_ARCHS4_coexpression', u'ZNF703_human_tf_ARCHS4_coexpression', u'ZFP91_human_tf_ARCHS4_coexpression', u'ZNF777_human_tf_ARCHS4_coexpression', u'ZBTB22_human_tf_ARCHS4_coexpression', u'ZBTB12_human_tf_ARCHS4_coexpression', u'TCF20_human_tf_ARCHS4_coexpression', u'TCF25_human_tf_ARCHS4_coexpression', u'ZC3H3_human_tf_ARCHS4_coexpression', u'NCOA3_human_tf_ARCHS4_coexpression', u'ZNF787_human_tf_ARCHS4_coexpression', u'ZMAT2_human_tf_ARCHS4_coexpression', u'FOXK1_human_tf_ARCHS4_coexpression', u'ZNF598_human_tf_ARCHS4_coexpression', u'ZNF768_human_tf_ARCHS4_coexpression', u'ZNF316_human_tf_ARCHS4_coexpression', u'PLXND1_human_tf_ARCHS4_coexpression', u'DVL2_human_tf_ARCHS4_coexpression', u'CSDE1_human_tf_ARCHS4_coexpression', u'AKAP8L_human_tf_ARCHS4_coexpression', u'BRD9_human_tf_ARCHS4_coexpression', u'HMG20B_human_tf_ARCHS4_coexpression', u'MBD3_human_tf_ARCHS4_coexpression', u'SRCAP_human_tf_ARCHS4_coexpression', u'ZC3H7B_human_tf_ARCHS4_coexpression', u'DOT1L_human_tf_ARCHS4_coexpression', u'SF3A2_human_tf_ARCHS4_coexpression', u'DVL1_human_tf_ARCHS4_coexpression', u'PLXNA1_human_tf_ARCHS4_coexpression', u'SLC2A4RG_human_tf_ARCHS4_coexpression', u'PRR12_human_tf_ARCHS4_coexpression', u'UBR4_human_tf_ARCHS4_coexpression', u'PLXNA3_human_tf_ARCHS4_coexpression', u'UNK_human_tf_ARCHS4_coexpression', u'SMAD6_human_tf_ARCHS4_coexpression', u'SMARCC2_human_tf_ARCHS4_coexpression', u'MINK1_human_tf_ARCHS4_coexpression', u'MBD1_human_tf_ARCHS4_coexpression', u'POGK_human_tf_ARCHS4_coexpression', u'SETDB1_human_tf_ARCHS4_coexpression', u'MTA2_human_tf_ARCHS4_coexpression', u'BRPF1_human_tf_ARCHS4_coexpression', u'ADAR_human_tf_ARCHS4_coexpression', u'DVL3_human_tf_ARCHS4_coexpression', u'PLXNB1_human_tf_ARCHS4_coexpression', u'ATMIN_human_tf_ARCHS4_coexpression', u'GATAD2A_human_tf_ARCHS4_coexpression', u'MECP2_human_tf_ARCHS4_coexpression', u'SSH3_human_tf_ARCHS4_coexpression', u'KAT5_human_tf_ARCHS4_coexpression', u'ANAPC2_human_tf_ARCHS4_coexpression', u'CCDC71_human_tf_ARCHS4_coexpression', u'TIGD5_human_tf_ARCHS4_coexpression', u'RBM10_human_tf_ARCHS4_coexpression', u'MTA1_human_tf_ARCHS4_coexpression', u'EP400_human_tf_ARCHS4_coexpression', u'H1FX_human_tf_ARCHS4_coexpression', u'H1F0_human_tf_ARCHS4_coexpression', u'NCOR2_human_tf_ARCHS4_coexpression', u'PLXNB2_human_tf_ARCHS4_coexpression', u'ZNF512B_human_tf_ARCHS4_coexpression', u'INF2_human_tf_ARCHS4_coexpression', u'REPIN1_human_tf_ARCHS4_coexpression'])]

	def testGenemap(self, gene, results):
		enrichr = Enrichr()
		ret = list(itertools.islice(enrichr.genemap(gene), 3))
		self.assertListEqual(ret, results)

	def dataProvider_testExport(self):
		yield GENES, 'Description', RESULTS_FILE.name, RESULTS

	def testExport(self, genes, description, retfile, results):
		self.maxDiff = None
		enrichr = Enrichr(cutoff = 1)
		enrichr.addList(genes, description)
		enrichr.enrich()
		enrichr.export(retfile)
		ret = []
		with open(retfile) as f:
			for line in f:
				line = line.strip()
				if not line or line.startswith('Term'): continue
				ret.append(Enrichr_Term._make(line.split('\t')))
		for i, result in enumerate(results):
			for field in result._fields:
				rval = getattr(result, field)
				oval = getattr(ret[i], field)
				if isinstance(rval, (int, float)):
					oval = float(oval)
					self.assertAlmostEqual(rval, oval, places=3)
				else:
					self.assertEqual(rval, oval)

	def dataProvider_testPlot(self):
		yield GENES, 'Description', path.join(tempfile.gettempdir(), 'plot.png')

	def testPlot(self, genes, description, plotfile):
		enrichr = Enrichr(cutoff = 1)
		enrichr.addList(genes, description)
		enrichr.enrich()
		enrichr.plot(plotfile)
		self.assertGreater(path.getsize(plotfile), 0)

if __name__ == '__main__':
	GENES_FILE.write('\n'.join(['%s\t%s' % (i, gene) for i, gene in enumerate(GENES)]).encode())
	GENES_FILE.close()
	testly.main(verbosity = 2)
