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

MT = 'Analyses/scz_allentries.mt'
PHENOTYPES_TABLE = 'data/scz_phenoinfo.txt'
IMPUTESEX_TABLE = 'Analyses/imputesex.ht'

SEXCHECK_LIST = 'Analyses/sexcheck.remove.sample_list'

PCA_LIST = 'Analyses/european.strict.sample_list'

IBD_SAMPLES = 'Analyses/ibd.remove.sample_list'
INITIAL_VARIANT_LIST = 'Analyses/initial_varlist_checkpoint.ht'
FINAL_VARIANT_LIST = 'Analyses/final_qc.keep.variant_list'

SAMPLE_BEFORE_QC_FILE = 'Analyses/final_qc.before.samples.tsv'
SAMPLE_AFTER_QC_FILE = 'Analyses/final_qc.after.samples.tsv'

sample_annotations = (hl.import_table(PHENOTYPES_TABLE).key_by('SAMPID'))
impute_sex_annotations = hl.read_table(IMPUTESEX_TABLE)

ht_pca_samples = hl.import_table(PCA_LIST, no_header=True, key='f0')
ht_ibd_samples = hl.import_table(IBD_SAMPLES, no_header=True, key='f0')
ht_sex_check_samples = hl.import_table(SEXCHECK_LIST, no_header=True, key='f0')

ht_initial_variants = hl.read_table(INITIAL_VARIANT_LIST)

ht_final_variants = hl.import_table(FINAL_VARIANT_LIST,
	types={'locus':hl.tlocus(reference_genome='GRCh37'), 'alleles':hl.tarray(hl.tstr)})

ht_final_variants = ht_final_variants.key_by(ht_final_variants.locus, ht_final_variants.alleles)

mt_before = hl.read_matrix_table(MT)
mt_before = mt_before.filter_cols(hl.is_defined(ht_pca_samples[mt_before.col_key]))

mt_before = mt_before.filter_cols(~hl.is_defined(ht_ibd_samples[mt_before.col_key]))
mt_before = mt_before.filter_cols(~hl.is_defined(ht_sex_check_samples[mt_before.col_key]))
mt_before = mt_before.filter_rows(hl.is_defined(ht_initial_variants[mt_before.row_key]))

mt_before = mt_before.annotate_cols(phenotype = sample_annotations[mt_before.col_key])
mt_before = mt_before.annotate_cols(imputesex = impute_sex_annotations[mt_before.col_key])

mt_before = hl.variant_qc(mt_before, name = 'qc')

mt_before = mt_before.annotate_rows(
	qc = mt_before.qc.annotate(AC=mt_before.qc.AC[1],
	AF = mt_before.qc.AF[1],
	homozygote_count = mt_before.qc.homozygote_count[1]))

mt_before = mt_before.filter_rows((mt_before.qc.AF > 0) & (mt_before.qc.AF < 1))
mt_before = hl.sample_qc(mt_before)

n = mt_before.count()

print('n samples:')
print(n[1])
print('n variants:')
print(n[0])

mt_before = mt_before.annotate_cols(sex = hl.case()
	.when(mt_before.imputesex.impute_sex.is_female, "Female")
	.default("Male"))

mt_after = mt_before.filter_rows(hl.is_defined(ht_final_variants[mt_before.row_key]))
mt_after = hl.sample_qc(mt_after)

n = mt_after.count()

print('n samples:')
print(n[1])
print('n variants:')
print(n[0])

mt_before.cols().select("sex", "phenotype", "sample_qc").flatten().export(SAMPLE_BEFORE_QC_FILE)
mt_after.cols().select("sex", "phenotype", "sample_qc").flatten().export(SAMPLE_AFTER_QC_FILE)
