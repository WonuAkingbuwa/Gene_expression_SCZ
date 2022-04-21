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

#export final plink file
MT_HARDCALLS = 'Analyses/scz_hardcalls.mt'
FINAL_PLINK_FILES = 'plink_prune/final_qc'

FINAL_SAMPLE_LIST = 'Analyses/final_qc.keep.sample_list'
FINAL_VARIANT_LIST = 'Analyses/final_qc.keep.variant_list'

ht_final_samples = hl.import_table(FINAL_SAMPLE_LIST, no_header=True, key='f0')
ht_final_variants = hl.import_table(FINAL_VARIANT_LIST, types={'locus':hl.tlocus(reference_genome='GRCh37'), 'alleles':hl.tarray(hl.tstr)})
ht_final_variants = ht_final_variants.key_by(ht_final_variants.locus, ht_final_variants.alleles)

mt = hl.read_matrix_table(MT_HARDCALLS)
mt = mt.filter_cols(hl.is_defined(ht_final_samples[mt.col_key]))

mt = mt.filter_rows(hl.is_defined(ht_final_variants[mt.row_key]))
mt = hl.variant_qc(mt, name='qc')
mt = mt.filter_rows((mt.qc.AF[0] > 0.05) & (mt.qc.AF[0] < 0.95) & ((mt.qc.call_rate > 0.98) | mt.locus.in_x_nonpar() | mt.locus.in_x_par())).persist()

for x in range(1,23):

        mt_chr = hl.filter_intervals(mt, [hl.parse_locus_interval(hl.eval(hl.str(x)), reference_genome='GRCh37')])
        n_chr = mt_chr.count_rows()

        print('\nn variants in chr')
        print(x)
        print(n_chr)

        hl.export_plink(mt_chr, FINAL_PLINK_FILES + '.chr' + str(x))

mt_chr = hl.filter_intervals(mt, [hl.parse_locus_interval('X', reference_genome='GRCh37')])
n_chr = mt_chr.count_rows()

print('\nn variants in chrX')
print(n_chr)

hl.export_plink(mt_chr, FINAL_PLINK_FILES + '.chr' + 'X')

#perform pruning in plink - go to plink script or run the command below
./15_prune_plink_final.sh
