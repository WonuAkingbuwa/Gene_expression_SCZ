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
import os 

####RV counts for each geneset####
QC_MT = 'Analyses/european.strict.mt'
RV_GENESETS = 'annotations/RVs_genesets/RVs_'
INTERVALS_GENESETS_FOLDER = 'annotations/geneset_bed_files/'

#read in mt and filter to autosomes 
mt = hl.read_matrix_table(QC_MT)
mt = mt.filter_rows(mt.locus.in_autosome())

#import intervals for genesets and write out resulting counts
resultfile = open("annotations/RVs_genesets/number_of_rv_variants_in_analyses.tsv", 'w')
for x in os.listdir(INTERVALS_GENESETS_FOLDER):
    if x.endswith('.bed'):
        INTERVALS_GENESETS_FILE = '{}/{}'.format(INTERVALS_GENESETS_FOLDER, x)
    else:
        continue
        
    intervals = hl.import_locus_intervals(INTERVALS_GENESETS_FILE, reference_genome='GRCh37', skip_invalid_intervals=True)
    mt_intervals = mt.annotate_rows(in_interval = hl.is_defined(intervals[mt.locus]))
    mt_intervals = mt_intervals.filter_rows(mt_intervals.in_interval)
    n_variants = mt_intervals.count_rows()

    print('\nn variants in filtered mt for')
    print(x)
    print(n_variants)
    resultfile.write("{}\t{}\n".format(x, n_variants))

    mt_intervals = mt_intervals.annotate_rows(is_singleton = hl.agg.sum(mt_intervals.GT.n_alt_alleles()) == 1) 
    mt_intervals = mt_intervals.annotate_rows(is_URV = (mt_intervals.is_singleton) & (mt_intervals.annotation.inGnomAD == True) & (~mt_intervals.annotation.indiscovEHR))
    mt_intervals = mt_intervals.filter_rows(~mt_intervals.is_URV) 

    mt_intervals = mt_intervals.filter_rows(mt_intervals.qc.AF[1] < 0.001) 

    mt_intervals = mt_intervals.annotate_cols(n_RV_SNP = hl.agg.count_where(mt_intervals.GT.is_non_ref() & hl.is_snp(mt_intervals.alleles[0], mt_intervals.alleles[1])),
                      n_RV_indel = hl.agg.count_where(mt_intervals.GT.is_non_ref() & hl.is_indel(mt_intervals.alleles[0], mt_intervals.alleles[1])),
                      n_coding_RV_SNP = hl.agg.count_where(mt_intervals.GT.is_non_ref() & hl.is_snp(mt_intervals.alleles[0], mt_intervals.alleles[1]) & (mt_intervals.annotation.consequence_category != "non_coding")),
                      n_coding_RV_indel = hl.agg.count_where(mt_intervals.GT.is_non_ref() & hl.is_indel(mt_intervals.alleles[0], mt_intervals.alleles[1]) & (mt_intervals.annotation.consequence_category != "non_coding")),
                      n_RV_PTV = hl.agg.count_where(mt_intervals.GT.is_non_ref() & (mt_intervals.annotation.consequence_category == "ptv")),
                      n_RV_damaging_missense = hl.agg.count_where(mt_intervals.GT.is_non_ref() & (mt_intervals.annotation.consequence_category == "damaging_missense")),
                      n_RV_other_missense = hl.agg.count_where(mt_intervals.GT.is_non_ref() & (mt_intervals.annotation.consequence_category == "other_missense")),
                      n_RV_synonymous = hl.agg.count_where(mt_intervals.GT.is_non_ref() & (mt_intervals.annotation.consequence_category == "synonymous")),
                      n_RV_non_coding = hl.agg.count_where(mt_intervals.GT.is_non_ref() & (mt_intervals.annotation.consequence_category == "non_coding")))

    mt_intervals.cols().flatten().export(RV_GENESETS + str(x) + '.tsv')

resultfile.close()

