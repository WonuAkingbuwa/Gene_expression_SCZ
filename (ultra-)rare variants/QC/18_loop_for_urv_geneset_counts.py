#cd to folder with scz mt data
cd /media/veracrypt10/Analyses/

#switch to hail 
conda activate hail

#run following command to try and prevent java heap space error - also loads ipython
PYSPARK_SUBMIT_ARGS="--driver-memory 16G --executor-memory 16G pyspark-shell" ipython

#to update hail  - run only when new versions are available 
#pip install hail --upgrade

#Import hail
import hail as hl
hl.init(min_block_size=128)

#import other python libraries that will be relevant later
from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()
import os 

####URV counts for each geneset####
QC_MT = 'Analyses/european.strict.mt'
URV_GENESETS = 'annotations/URVs_genesets/URVs_'
INTERVALS_GENESETS_FOLDER = 'annotations/geneset_bed_files/'

#read in mt and filter to autosomes 
mt = hl.read_matrix_table(QC_MT)
mt = mt.filter_rows(mt.locus.in_autosome())

#import intervals for genesets and write out resulting counts
resultfile = open("annotations/URVs_genesets/number_of_urv_variants_in_analyses.tsv", 'w')
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

    mt_intervals = mt_intervals.filter_rows((mt_intervals.is_singleton) & (mt_intervals.annotation.inGnomAD == True) & (~mt_intervals.annotation.indiscovEHR)) 

    mt_intervals = mt_intervals.annotate_cols(n_URV_SNP = hl.agg.count_where(mt_intervals.GT.is_non_ref() & hl.is_snp(mt_intervals.alleles[0], mt_intervals.alleles[1])),
                      n_URV_indel = hl.agg.count_where(mt_intervals.GT.is_non_ref() & hl.is_indel(mt_intervals.alleles[0], mt_intervals.alleles[1])),
                      n_coding_URV_SNP = hl.agg.count_where(mt_intervals.GT.is_non_ref() & hl.is_snp(mt_intervals.alleles[0], mt_intervals.alleles[1]) & (mt_intervals.annotation.consequence_category != "non_coding")),
                      n_coding_URV_indel = hl.agg.count_where(mt_intervals.GT.is_non_ref() & hl.is_indel(mt_intervals.alleles[0], mt_intervals.alleles[1]) & (mt_intervals.annotation.consequence_category != "non_coding")),
                      n_URV_PTV = hl.agg.count_where(mt_intervals.GT.is_non_ref() & (mt_intervals.annotation.consequence_category == "ptv")),
                      n_URV_damaging_missense = hl.agg.count_where(mt_intervals.GT.is_non_ref() & (mt_intervals.annotation.consequence_category == "damaging_missense")),
                      n_URV_other_missense = hl.agg.count_where(mt_intervals.GT.is_non_ref() & (mt_intervals.annotation.consequence_category == "other_missense")),
                      n_URV_synonymous = hl.agg.count_where(mt_intervals.GT.is_non_ref() & (mt_intervals.annotation.consequence_category == "synonymous")),
                      n_URV_non_coding = hl.agg.count_where(mt_intervals.GT.is_non_ref() & (mt_intervals.annotation.consequence_category == "non_coding")))

    mt_intervals.cols().flatten().export(URV_GENESETS + str(x) + '.tsv')

resultfile.close()
