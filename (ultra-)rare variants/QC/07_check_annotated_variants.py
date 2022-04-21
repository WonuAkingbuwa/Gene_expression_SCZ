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

####check number of relevant annotated variants####
ANNOTATION_TABLE = 'annotations/gene.ht'

ht = hl.read_table(ANNOTATION_TABLE)
ht = ht.annotate(type = hl.case().when((ht.alleles[0].length() == 1) | (ht.alleles[1].length() == 1), "SNP")
	.when(ht.alleles[0].length() < ht.alleles[1].length(), "Insertion")
	.when(ht.alleles[0].length() > ht.alleles[1].length(), "Deletion")
	.default("No type"))

n = ht.count()
n_inGnomAD = ht.aggregate(hl.agg.counter(ht.inGnomAD))
n_indiscovEHR = ht.aggregate(hl.agg.counter(ht.indiscovEHR))
n_have_gene = ht.aggregate(hl.agg.counter(hl.is_defined(ht.vep.worst_csq_for_variant_canonical.gene_symbol)))
n_consequence_category = ht.aggregate(hl.agg.counter(ht.consequence_category))
n_most_severe_consequence = ht.aggregate(hl.agg.counter(ht.vep.worst_csq_for_variant_canonical.most_severe_consequence))

print('n variants:')
print(n) 

print('gnomAD')
print(n_inGnomAD)

print('discovEHR')
print(n_indiscovEHR) 

print('have_gene')
print(n_have_gene) 

print('consequence_category')
print(n_consequence_category) 

print('most_severe_consequence')
print(n_most_severe_consequence) 

# Note that this is not restricted to the prefiltered variants.
# Let's check the counts filtered to the target intervals:
ht_initial_variants = hl.import_table('Analyses/prefilter.keep.variant_list', types={'locus':hl.tlocus(reference_genome='GRCh37'), 'alleles':hl.tarray(hl.tstr)})
ht_initial_variants = ht_initial_variants.key_by(locus=ht_initial_variants.locus, alleles=ht_initial_variants.alleles)
ht = ht.filter(hl.is_defined(ht_initial_variants[ht.key]))

n = ht.count()
n_inGnomAD = ht.aggregate(hl.agg.counter(ht.inGnomAD))
n_indiscovEHR = ht.aggregate(hl.agg.counter(ht.indiscovEHR))
n_have_gene = ht.aggregate(hl.agg.counter(hl.is_defined(ht.vep.worst_csq_for_variant_canonical.gene_symbol)))
n_consequence_category = ht.aggregate(hl.agg.counter(ht.consequence_category))
n_most_severe_consequence = ht.aggregate(hl.agg.counter(ht.vep.worst_csq_for_variant_canonical.most_severe_consequence))

print('n variants:')
print(n) 

print('gnomAD')
print(n_inGnomAD) 

print('discovEHR')
print(n_indiscovEHR) 

print('have_gene')
print(n_have_gene) 

print('consequence_category')
print(n_consequence_category) 

print('most_severe_consequence')
print(n_most_severe_consequence)
