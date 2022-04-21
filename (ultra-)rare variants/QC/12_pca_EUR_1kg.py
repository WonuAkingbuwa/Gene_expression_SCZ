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
MT_1KG = 'Analyses/1KG_autosomes_GRCh37.mt'
PCA_SCORES_EUR = 'Analyses/pca_scores.european.1kg.tsv'

PHENOTYPES_TABLE = 'data/scz_phenoinfo.txt'
POPULATIONS_1KG = 'data/1kg_annotations.txt'
INITIAL_SAMPLES = 'Analyses/initial_qc.keep.sample_list'
PRUNED_VARIANTS = 'plink_prune/pruned/scz_prune.keep.variant_list'
IBD_SAMPLES = 'Analyses/ibd.remove.sample_list'

EUROPEAN_SAMPLES = 'Analyses/european.strict.sample_list'

mt_1kg = hl.read_matrix_table(MT_1KG)
mt_1kg = hl.split_multi_hts(mt_1kg)
# This is to enable a join later.
mt_1kg = mt_1kg.select_entries("GT")

populations_1kg = (hl.import_table(POPULATIONS_1KG).key_by('Sample'))
mt_1kg = mt_1kg.annotate_cols(pop_1kg = populations_1kg[mt_1kg.s])
mt_1kg = mt_1kg.filter_cols(mt_1kg.pop_1kg.SuperPopulation == "EUR")
mt_1kg = mt_1kg.select_cols()

ht_initial_samples = hl.import_table(INITIAL_SAMPLES, no_header=True, key='f0')
ht_pruned_variants = hl.import_table(PRUNED_VARIANTS, no_header=True)

ht_pruned_variants = ht_pruned_variants.annotate(**hl.parse_variant(ht_pruned_variants.f0, reference_genome='GRCh37'))
ht_pruned_variants = ht_pruned_variants.key_by(ht_pruned_variants.locus, ht_pruned_variants.alleles)

ht_ibd_samples = hl.import_table(IBD_SAMPLES, no_header=True, key='f0')
ht_european_samples = hl.import_table(EUROPEAN_SAMPLES, no_header=True, key='f0')

mt = hl.read_matrix_table(MT_HARDCALLS)

mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.filter_cols(~hl.is_defined(ht_ibd_samples[mt.col_key]))
mt = mt.filter_cols(hl.is_defined(ht_european_samples[mt.col_key]))
mt = mt.filter_rows(hl.is_defined(ht_pruned_variants[mt.row_key]))

mt = mt.union_cols(mt_1kg)

n = mt.count()

print('n samples:')
print(n[1])
print('n variants:')
print(n[0])

pca_output = hl.hwe_normalized_pca(mt.GT, k=10)
pca_output = pca_output[1].key_by('s')

sample_annotations = (hl.import_table(PHENOTYPES_TABLE).key_by('SAMPID'))

pca_output = pca_output.annotate(phenotype = sample_annotations[pca_output.s])
pca_output = pca_output.annotate(
	PC1 = pca_output.scores[0],
	PC2 = pca_output.scores[1],
	PC3 = pca_output.scores[2],
	PC4 = pca_output.scores[3],
	PC5 = pca_output.scores[4],
	PC6 = pca_output.scores[5],
	PC7 = pca_output.scores[6],
	PC8 = pca_output.scores[7],
	PC9 = pca_output.scores[8],
	PC10 = pca_output.scores[9])

pca_output = pca_output.annotate(pop_1kg = populations_1kg[pca_output.s])

pca_output = pca_output.annotate(ANALYSIS_CAT = hl.case().when(hl.is_defined(pca_output.phenotype.ANALYSIS_CAT), pca_output.phenotype.ANALYSIS_CAT)
                                 .default("1KG"),
                                 SITE = hl.case().when(hl.is_defined(pca_output.phenotype.SITE), pca_output.phenotype.SITE)
                                 .default("1KG"),
                                 SAMP_SOURCE = hl.case().when(hl.is_defined(pca_output.phenotype.SAMP_SOURCE), pca_output.phenotype.SAMP_SOURCE)
                                 .default("1KG"),
                                 PRIMARY_DISEASE = hl.case().when(hl.is_defined(pca_output.phenotype.PRIMARY_DISEASE), pca_output.phenotype.PRIMARY_DISEASE)
                                 .default("1KG"),
                                 SUPER_POPULATION = hl.case().when(hl.is_defined(pca_output.pop_1kg.SuperPopulation), pca_output.pop_1kg.SuperPopulation)
                                 .default(pca_output.phenotype.ANALYSIS_CAT),
                                 POPULATION = hl.case().when(hl.is_defined(pca_output.pop_1kg.Population), pca_output.pop_1kg.Population)
                                 .default(pca_output.phenotype.ANALYSIS_CAT)
).repartition(128).persist()

pca_output.flatten().export(output=PCA_SCORES_EUR)
