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

MT_HARDCALLS = 'Analyses/scz_hardcalls.mt'
PCA_SCORES = 'Analyses/pca_scores.tsv'

PHENOTYPES_TABLE = 'data/scz_phenoinfo.txt'
INITIAL_SAMPLES = 'Analyses/initial_qc.keep.sample_list'
PRUNED_VARIANTS = 'plink_prune/pruned/scz_prune.keep.variant_list'
IBD_SAMPLES = 'Analyses/ibd.remove.sample_list'

mt = hl.read_matrix_table(MT_HARDCALLS)
sample_annotations = (hl.import_table(PHENOTYPES_TABLE).key_by('SAMPID'))

ht_initial_samples = hl.import_table(INITIAL_SAMPLES, no_header=True, key='f0')
ht_pruned_variants = hl.import_table(PRUNED_VARIANTS, no_header=True)
ht_ibd_samples = hl.import_table(IBD_SAMPLES, no_header=True, key='f0')

ht_pruned_variants = ht_pruned_variants.annotate(**hl.parse_variant(ht_pruned_variants.f0, reference_genome='GRCh37'))
ht_pruned_variants = ht_pruned_variants.key_by(ht_pruned_variants.locus, ht_pruned_variants.alleles)

mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.filter_cols(~hl.is_defined(ht_ibd_samples[mt.col_key]))
mt = mt.filter_rows(hl.is_defined(ht_pruned_variants[mt.row_key]))

mt = mt.annotate_cols(phenotype = sample_annotations[mt.s]).repartition(128).persist()

n = mt.count()

print('n samples:')
print(n[1])
print('n variants:')
print(n[0])

pca_output = hl.hwe_normalized_pca(mt.GT, k=10)
pca_output = pca_output[1].key_by('s')

pca_output = pca_output.annotate(phenotype = sample_annotations[pca_output.s]).repartition(128).persist()

pca_output = pca_output.annotate(PC1 = pca_output.scores[0],
	PC2 = pca_output.scores[1], PC3 = pca_output.scores[2],
	PC4 = pca_output.scores[3], PC5 = pca_output.scores[4],
	PC6 = pca_output.scores[5], PC7 = pca_output.scores[6],
	PC8 = pca_output.scores[7], PC9 = pca_output.scores[8],
	PC10 = pca_output.scores[9])

pca_output.flatten().export(output=PCA_SCORES)
