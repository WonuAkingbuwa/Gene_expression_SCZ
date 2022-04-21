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

#assign output file names
PLINK_FILES = 'plink_prune/filterGT'

#import initial sample qc list
ht_initial_samples = hl.import_table('Analyses/initial_qc.keep.sample_list', no_header=True, key='f0')

#import initial variant list
ht_initial_variants = hl.import_table('Analyses/prefilter.keep.variant_list', types={'locus':hl.tlocus(reference_genome='GRCh37'), 'alleles':hl.tarray(hl.tstr)})
ht_initial_variants = ht_initial_variants.key_by(ht_initial_variants.locus, ht_initial_variants.alleles)

#import hardcalls mt
mt = hl.read_matrix_table('Analyses/scz_hardcalls.mt')
mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))

mt = mt.filter_rows(hl.is_defined(ht_initial_variants[mt.row_key]))
mt = mt.filter_rows(mt.locus.in_x_nonpar() | mt.locus.in_autosome_or_par())
mt = hl.variant_qc(mt, name = 'variant_qc')
mt = mt.filter_rows((mt.variant_qc.AF[0] > 0.05) & (mt.variant_qc.AF[0] < 0.95) & ((mt.variant_qc.call_rate > 0.98) | mt.locus.in_x_nonpar() | mt.locus.in_x_par())).persist()

n = mt.count()

pprint('n samples:')
print(n[1]) 
pprint('n variants:')
print(n[0]) 

#write out dataset for performing pruning, in hail mt format
#mt.repartition(100, shuffle=False).write('Analyses/scz_qc_prune.mt', overwrite=True)

#read in mt that will be used to make plink files for pruning - useful if rerunning later, no need to go through previous steps
#mt = hl.read_matrix_table('Analyses/scz_qc_prune.mt')

#write file for plink - chr(x) format appears to be invalid for GRCh37 so switching to just numbers. use hashed out version if reference genome is build 38
for x in range(1,23):

        mt_chr = hl.filter_intervals(mt, [hl.parse_locus_interval(hl.eval(hl.str(x)), reference_genome='GRCh37')])
        n_chr = mt_chr.count_rows()
        print('\nn variants in chr')
        print(x)
        print(n_chr)
        hl.export_plink(mt_chr, PLINK_FILES + '.chr' + str(x))

mt_chr = hl.filter_intervals(mt, [hl.parse_locus_interval('X', reference_genome='GRCh37')])
n_chr = mt_chr.count_rows()

print('\nn variants in chrX')
print(n_chr)

hl.export_plink(mt_chr, PLINK_FILES + '.chr' + 'X')

####perform pruning in plink - go to plink script or run the command below####
./04_prune_plink.sh
