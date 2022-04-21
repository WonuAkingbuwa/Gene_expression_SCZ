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

####URV counts####
QC_MT = 'Analyses/european.strict.mt'
RV_FILE = 'annotations/RVs_after_QC.tsv'

#read in mt and filter to autosomes 
mt = hl.read_matrix_table(QC_MT)
mt = mt.filter_rows(mt.locus.in_autosome())

mt = mt.annotate_rows(is_singleton = hl.agg.sum(mt.GT.n_alt_alleles()) == 1)
mt = mt.annotate_rows(is_URV = (mt.is_singleton) & (mt.annotation.inGnomAD == True) & (~mt.annotation.indiscovEHR))

#count number of URVs and singletons
n = mt.count()

print('')
print('n samples:')
print(n[1]) 
print('n variants:')
print(n[0]) 

n_issingleton = mt.aggregate_rows(hl.agg.counter(mt.is_singleton))
print('is_singleton')
print(n_issingleton) 

n_isURV = mt.aggregate_rows(hl.agg.counter(mt.is_URV))
print('is_URV')
print(n_isURV) 

#filter out URVs
mt = mt.filter_rows(~mt.is_URV)
mt.is_URV.show() 

n = mt.count()

print('')
print('n samples:')
print(n[1]) 
print('n variants:')
print(n[0]) 

mt = mt.filter_rows(mt.qc.AF[1] < 0.001) 

n = mt.count()

print('')
print('n samples:')
print(n[1]) 
print('n variants:')
print(n[0]) 

mt = mt.annotate_cols(n_RV_SNP = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_snp(mt.alleles[0], mt.alleles[1])),
                      n_RV_indel = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_indel(mt.alleles[0], mt.alleles[1])),
                      n_coding_RV_SNP = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_snp(mt.alleles[0], mt.alleles[1]) & (mt.annotation.consequence_category != "non_coding")),
                      n_coding_RV_indel = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_indel(mt.alleles[0], mt.alleles[1]) & (mt.annotation.consequence_category != "non_coding")),
                      n_RV_PTV = hl.agg.count_where(mt.GT.is_non_ref() & (mt.annotation.consequence_category == "ptv")),
                      n_RV_damaging_missense = hl.agg.count_where(mt.GT.is_non_ref() & (mt.annotation.consequence_category == "damaging_missense")),
                      n_RV_other_missense = hl.agg.count_where(mt.GT.is_non_ref() & (mt.annotation.consequence_category == "other_missense")),
                      n_RV_synonymous = hl.agg.count_where(mt.GT.is_non_ref() & (mt.annotation.consequence_category == "synonymous")),
                      n_RV_non_coding = hl.agg.count_where(mt.GT.is_non_ref() & (mt.annotation.consequence_category == "non_coding")))

mt.cols().flatten().export(RV_FILE)
