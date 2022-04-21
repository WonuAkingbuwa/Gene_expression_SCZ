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
IBD_OUTPUT = 'Analyses/ibd.tsv'

PHENOTYPES_TABLE = 'data/scz_phenoinfo.txt'
PRUNED_VARIANTS = 'plink_prune/pruned/scz_prune.keep.variant_list'

INITIAL_SAMPLES = 'Analyses/initial_qc.keep.sample_list'

ht_initial_samples = hl.import_table(INITIAL_SAMPLES, no_header=True, key='f0')
ht_pruned_variants = hl.import_table(PRUNED_VARIANTS, no_header=True)

ht_pruned_variants = ht_pruned_variants.annotate(**hl.parse_variant(ht_pruned_variants.f0, reference_genome='GRCh37'))
ht_pruned_variants = ht_pruned_variants.key_by(ht_pruned_variants.locus, ht_pruned_variants.alleles)

sample_annotations = (hl.import_table(PHENOTYPES_TABLE).key_by('SAMPID'))

mt = hl.read_matrix_table(MT_HARDCALLS)
mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.filter_rows(hl.is_defined(ht_pruned_variants[mt.row_key]))
mt = mt.annotate_cols(phenotype = sample_annotations[mt.s]).repartition(128).persist()

n = mt.count()

print('n samples:')
print(n[1])
print('n variants:')
print(n[0])

ibd_table = hl.identity_by_descent(mt, min=0.1)
ibd_table = ibd_table.flatten()
ibd_table.export(IBD_OUTPUT)

#move to ibd_filer/ibd_plot.R script 
