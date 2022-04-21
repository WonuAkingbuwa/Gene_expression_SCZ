#switch to hail 
conda activate hail
PYSPARK_SUBMIT_ARGS="--driver-memory 16G --executor-memory 16G pyspark-shell" ipython

#Import hail
import hail as hl
hl.init(min_block_size=128)

#import other python libraries that will be relevant later
from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

QC_MT = 'Analyses/european.strict.mt'
QC_HARDCALLS_MT = 'Analyses/european.strict.hardcalls.mt'
PHENOTYPES_TABLE = 'data/scz_phenoinfo.txt'

FINAL_VARIANT_QC_FILE = 'Analyses/final_qc.variants.tsv.gz'
FINAL_SAMPLE_QC_FILE = 'Analyses/final_qc.samples.tsv.gz'

mt = hl.read_matrix_table(QC_HARDCALLS_MT)
ht = mt.cols()

ht.describe() #use this to determine what can be selected in the subsequent sections 

ht.select(
	SAMP_SOURCE = ht.phenotype.SAMP_SOURCE,
	ANALYSIS_CAT = ht.phenotype.ANALYSIS_CAT,
	PRIMARY_DISEASE = ht.phenotype.PRIMARY_DISEASE,
	IS_FEMALE = ht.imputesex.impute_sex.is_female,
	SEX = ht.phenotype.SEX,
	QC = ht.sample_qc,
	PCA = ht.pca).flatten().export(FINAL_SAMPLE_QC_FILE)

# Want to write out the alternative allele count!
ht = mt.rows()

ht.describe() #use this to determine what can be selected in the subsequent sections

ht.select(
	rsid = ht.annotation.rsid,
	was_split = ht.annotation.was_split,
	gene_symbol = ht.annotation.vep.worst_csq_for_variant_canonical.gene_symbol,
	gene_id = ht.annotation.vep.worst_csq_for_variant_canonical.gene_id,

	# can add other annotations later on should they be necessary for any reason

	most_severe_consequence = ht.annotation.vep.worst_csq_for_variant_canonical.most_severe_consequence,
	consequence_category = ht.annotation.consequence_category,

	sift_prediction = ht.annotation.vep.worst_csq_for_variant_canonical.sift_pred,
	polyphen2_hdiv_prediction = ht.annotation.vep.worst_csq_for_variant_canonical.polyphen2_hdiv_pred,
	polyphen2_hvar_prediction = ht.annotation.vep.worst_csq_for_variant_canonical.polyphen2_hvar_pred,
	lrt_prediction = ht.annotation.vep.worst_csq_for_variant_canonical.lrt_pred,
	mutationtaster_prediction = ht.annotation.vep.worst_csq_for_variant_canonical.mutationtaster_pred,
	mutationassessor_prediction = ht.annotation.vep.worst_csq_for_variant_canonical.mutationassessor_pred,
	provean_prediction = ht.annotation.vep.worst_csq_for_variant_canonical.provean_pred,


	inGnomAD = ht.annotation.inGnomAD,
	indiscovEHR = ht.annotation.indiscovEHR,
	mpc_score = ht.annotation.mpc.MPC,

	qc = ht.qc).flatten().export(FINAL_VARIANT_QC_FILE)
