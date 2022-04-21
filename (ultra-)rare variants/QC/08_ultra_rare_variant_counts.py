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
MT = 'Analyses/scz_allentries.mt'
ANNOTATION_TABLE = 'annotations/gene.ht'
URV_FILE = 'annotations/URVs.tsv'
PHENOTYPES_TABLE = 'data/scz_phenoinfo.txt'

INITIAL_SAMPLES = 'Analyses/initial_qc.keep.sample_list'
INITIAL_VARIANT_LIST = 'Analyses/prefilter.keep.variant_list'
SEXCHECK_SAMPLES = 'Analyses/sexcheck.remove.sample_list'
IBD_SAMPLES = 'Analyses/ibd.remove.sample_list'

ht_initial_variants = hl.import_table(INITIAL_VARIANT_LIST,
    types={'locus':hl.tlocus(reference_genome='GRCh37'), 'alleles':hl.tarray(hl.tstr)})
ht_initial_variants = ht_initial_variants.key_by(locus=ht_initial_variants.locus, alleles=ht_initial_variants.alleles)
ht_initial_samples = hl.import_table(INITIAL_SAMPLES, no_header=True, key='f0')
ht_ibd_samples = hl.import_table(IBD_SAMPLES, no_header=True, key='f0')
ht_sexcheck_samples = hl.import_table(SEXCHECK_SAMPLES, no_header=True, key='f0')

mt = hl.read_matrix_table(MT)
mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))  
mt = mt.filter_cols(~hl.is_defined(ht_ibd_samples[mt.col_key]))
mt = mt.filter_cols(~hl.is_defined(ht_sexcheck_samples[mt.col_key]))
mt = mt.filter_rows(hl.is_defined(ht_initial_variants[mt.row_key]))
# Drop some fields that are not needed.
mt = mt.drop('a_index', 'qual', 'rsid', 'info', 'filters', 'was_split')

sample_annotations = (hl.import_table(PHENOTYPES_TABLE).key_by('SAMPID'))
constraint_annotations = hl.read_table(ANNOTATION_TABLE)

mt = mt.annotate_cols(phenotype = sample_annotations[mt.s])
mt = mt.annotate_rows(constraint = constraint_annotations[mt.row_key])
mt = mt.annotate_rows(is_singleton = hl.agg.sum(mt.GT.n_alt_alleles()) == 1)

# Filter variants based on Gnomad/discovEHR variant lists.
mt = mt.filter_rows((mt.is_singleton) & (mt.constraint.inGnomAD == True) & (~mt.constraint.indiscovEHR))

# Annotate variants with LoF/damaging missense annotation.
mt = mt.annotate_cols(n_URV_SNP = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_snp(mt.alleles[0], mt.alleles[1])),
                      n_URV_indel = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_indel(mt.alleles[0], mt.alleles[1])),
                      n_coding_URV_SNP = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_snp(mt.alleles[0], mt.alleles[1]) & (mt.constraint.consequence_category != "non_coding")),
                      n_coding_URV_indel = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_indel(mt.alleles[0], mt.alleles[1]) & (mt.constraint.consequence_category != "non_coding")),
                      n_URV_PTV = hl.agg.count_where(mt.GT.is_non_ref() & (mt.constraint.consequence_category == "ptv")),
                      n_URV_damaging_missense = hl.agg.count_where(mt.GT.is_non_ref() & (mt.constraint.consequence_category == "damaging_missense")),
                      n_URV_other_missense = hl.agg.count_where(mt.GT.is_non_ref() & (mt.constraint.consequence_category == "other_missense")),
                      n_URV_synonymous = hl.agg.count_where(mt.GT.is_non_ref() & (mt.constraint.consequence_category == "synonymous")),
                      n_URV_non_coding = hl.agg.count_where(mt.GT.is_non_ref() & (mt.constraint.consequence_category == "non_coding")))

mt.cols().flatten().export(URV_FILE)


