#switch to hail 
conda activate hail
PYSPARK_SUBMIT_ARGS="--driver-memory 16G --executor-memory 16G pyspark-shell" ipython

####Import hail####
import hail as hl
hl.init(min_block_size=128)

#import other python libraries that will be relevant later
from hail.plot import show
from pprint import pprint
hl.plot.output_notebook() 

MT_HARDCALLS = 'Analyses/scz_hardcalls.mt' 

IMPUTESEX_TABLE = 'Analyses/imputesex.ht'
IMPUTESEX_FILE = 'Analyses/imputesex.tsv'
Y_NCALLED = 'Analyses/ycalled.tsv'

INITIAL_SAMPLES = 'Analyses/initial_qc.keep.sample_list'
PRUNED_CHRX_VARIANTS = 'plink_prune/pruned/scz_chrX.prune.in'

PHENOTYPES_TABLE = 'data/scz_phenoinfo.txt'

ht_initial_samples = hl.import_table(INITIAL_SAMPLES, no_header=True, key='f0')
ht_pruned_chrx_variants = hl.import_table(PRUNED_CHRX_VARIANTS, no_header=True)
sample_annotations = (hl.import_table(PHENOTYPES_TABLE).key_by('SAMPID'))

ht_pruned_chrx_variants = ht_pruned_chrx_variants.annotate(**hl.parse_variant(ht_pruned_chrx_variants.f0, reference_genome='GRCh37'))
ht_pruned_chrx_variants = ht_pruned_chrx_variants.key_by(ht_pruned_chrx_variants.locus, ht_pruned_chrx_variants.alleles)

mt = hl.read_matrix_table(MT_HARDCALLS)
mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.filter_rows(hl.is_defined(ht_pruned_chrx_variants[mt.row_key]))

n = mt.count()

print('n samples:')
print(n[1]) 
print('n variants:')
print(n[0])

imputed_sex = hl.impute_sex(mt.GT, female_threshold=0.6, male_threshold=0.6)
mt = mt.annotate_cols(phenotype = sample_annotations[mt.s])
mt = mt.annotate_cols(impute_sex = imputed_sex[mt.s])

mt.cols().select('impute_sex', 'phenotype').flatten().export(IMPUTESEX_FILE)
mt.cols().write(IMPUTESEX_TABLE, overwrite=True)

# Determine non-missing allele count on the y.
mt = hl.read_matrix_table(MT_HARDCALLS)
mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.filter_rows(mt.locus.in_y_nonpar() | mt.locus.in_y_par())
mt = hl.sample_qc(mt, name='qc')

mt_cols = mt.cols()
mt_cols.select(n_called=mt_cols.qc.n_called).export(Y_NCALLED)

#move to impute_sex.R script 
