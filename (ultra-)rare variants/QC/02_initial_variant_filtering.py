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

#read in hardcalls file 
mt = hl.read_matrix_table('Analyses/scz_hardcalls.mt')

#remove variants that fail VQSR, i.e. select only "PASS" variants (others are likely false positives) - those that are blank in the filter field (https://hail.is/docs/0.2/methods/impex.html#hail.methods.import_vcf)
#annotate variants with flag indicating if they have failed VQSR.
mt = mt.annotate_rows(fail_VQSR = hl.len(mt.filters) != 0)
fail_VQSR = mt.filter_rows(mt.fail_VQSR).count_rows()
print('n variants failing VQSR:')
pprint(fail_VQSR)

#export variant annotations and variant keytable
mt_rows = mt.rows()
mt_rows.select(mt_rows.fail_VQSR).export('Analyses/prefilter_metrics.tsv')

mt = mt.filter_rows(mt.fail_VQSR, keep = False)

intervals = [hl.parse_locus_interval(x, reference_genome='GRCh37') for x in ['1:START-22:END', 'X:START-X:END', 'Y:START-Y:END']]
mt = hl.filter_intervals(mt, intervals)

#filter out invariant rows
mt = hl.variant_qc(mt)
mt = mt.filter_rows((mt.variant_qc.AF[0] > 0.0) & (mt.variant_qc.AF[0] < 1.0))

mt_rows_filter = mt.rows().select().export('Analyses/prefilter.keep.variant_list')

n_variants = hl.import_table('Analyses/prefilter.keep.variant_list').count()

print('n variants after initial filter:')
print(n_variants)

#after this you can run initial_variant_filter_info.R
