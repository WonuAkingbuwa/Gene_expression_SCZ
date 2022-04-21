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

#import filtered entries files
mt = hl.read_matrix_table('Analyses/scz_allentries.mt')

#import variant keep list
variants_to_filter = hl.import_table('Analyses/prefilter.keep.variant_list',
	types={'locus':hl.tlocus(reference_genome='GRCh37'), 'alleles':hl.tarray(hl.tstr)})
variants_to_filter = variants_to_filter.key_by(locus=variants_to_filter.locus, alleles=variants_to_filter.alleles)

#import phenotype table and annotate to genotype files - includes those with bipolar
phenotypes = hl.import_table('data/scz_phenoinfo.txt')
sample_annotations = phenotypes.key_by('SAMPID') #will aid annotation to the scz mt in subsequent steps

#filter data to variants we want to keep and annotate to phenotype file
mt = mt.filter_rows(hl.is_defined(variants_to_filter[mt.row_key]))
mt = mt.annotate_cols(phenotype = sample_annotations[mt.s])

n = mt.count()

pprint('n samples:')
print(n[1])
pprint('n variants:')
print(n[0])

mt = hl.sample_qc(mt, name = 'sample_qc')
mt.cols().select('phenotype', 'sample_qc').flatten().export(output='Analyses/initial_sample_qc.tsv') #including 'sample_qc' in the select command allows you to also export sample qc metrics so that they can be checked in R. Easier to manipulate than doing it in hail

##followed by initial_sample_qc_filter.R and initial_sample_qc_plots.R for further sample qc


