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

#import all entries filtered file 
MT = 'Analyses/scz_allentries.mt'
GNOMAD_SITES_MT = 'annotations/gnomad.exomes.r2.1.1.sites.mt'
DISCOVHER_VARIANTS = 'annotations/discovEHR.mt'

mt = hl.read_matrix_table(MT)

#import initial sample qc list
ht_initial_samples = hl.import_table('Analyses/initial_qc.keep.sample_list', no_header=True, key='f0')

#import initial variant list
ht_initial_variants = hl.import_table('Analyses/prefilter.keep.variant_list', types={'locus':hl.tlocus(reference_genome='GRCh37'), 'alleles':hl.tarray(hl.tstr)})
ht_initial_variants = ht_initial_variants.key_by(ht_initial_variants.locus, ht_initial_variants.alleles)

#filter by variant and sample keep lists
mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.filter_rows(hl.is_defined(ht_initial_variants[mt.row_key]))

#filter on rare variants 
mt = hl.variant_qc(mt)
mt.describe()
mt = mt.filter_rows(mt.variant_qc.AF[1] < 0.001)

mt.rows().select().export('Analyses/prefilter.keep.rare_variant_list')

n_variants = hl.import_table('Analyses/prefilter.keep.rare_variant_list').count()

print('n variants after initial rare variant filter:')
print(n_variants)

#discovehr
discovehr_variants_mt = hl.read_matrix_table(DISCOVHER_VARIANTS)
discovehr_variants_mt.describe()
discovehr_variants_mt.info.AF_LT_001.show()
discovehr_rare_variants_mt = discovehr_variants_mt.filter_rows(discovehr_variants_mt.info.AF_LT_001 == True)
discovehr_rare_variants_mt.show() 
discovehr_rare_variants_mt.count_rows()
discovehr_rare_variants_mt.rows().select().export('Analyses/discovehr.keep.rare_variant_list')

#gnomad - change all references to non-neuro data here
gnomad_mt = hl.read_matrix_table(GNOMAD_SITES_MT)
gnomad_mt.describe()
gnomad_mt.info.AF.show()
gnomad_rare_variants_mt = gnomad_mt.filter_rows(gnomad_mt.info.AF[0] < 0.001)
gnomad_rare_variants_mt.info.AF.show()
gnomad_rare_variants_mt.count_rows() 
gnomad_rare_variants_mt.rows().select().export('Analyses/gnomad.keep.rare_variant_mt_list')

#compare with ht list

####RV counts - rare variants defined using just scz dataset i.e. including URVs####
MT = 'Analyses/scz_allentries.mt'
ANNOTATION_TABLE = 'annotations/gene.ht'
RV_FILE = 'annotations/RVs.tsv'
PHENOTYPES_TABLE = 'data/scz_phenoinfo.txt'

INITIAL_SAMPLES = 'Analyses/initial_qc.keep.sample_list'
INITIAL_RARE_VARIANT_LIST = 'Analyses/prefilter.keep.rare_variant_list'
#GNOMAD_RARE_VARIANT_LIST = 'Analyses/gnomad.keep.rare_variant_list'
#DISCOVEHR_RARE_VARIANT_LIST = 'Analyses/discovehr.keep.rare_variant_list'
SEXCHECK_SAMPLES = 'Analyses/sexcheck.remove.sample_list'
IBD_SAMPLES = 'Analyses/ibd.remove.sample_list'

ht_initial_rare_variants = hl.import_table(INITIAL_RARE_VARIANT_LIST, types={'locus':hl.tlocus(reference_genome='GRCh37'), 'alleles':hl.tarray(hl.tstr)})
ht_initial_rare_variants = ht_initial_rare_variants.key_by(ht_initial_rare_variants.locus, ht_initial_rare_variants.alleles)

#ht_gnomad_rare_variants = hl.import_table(GNOMAD_RARE_VARIANT_LIST, types={'locus':hl.tlocus(reference_genome='GRCh37'), 'alleles':hl.tarray(hl.tstr)})
#ht_gnomad_rare_variants = ht_gnomad_rare_variants.key_by(ht_gnomad_rare_variants.locus, ht_gnomad_rare_variants.alleles)

#ht_discovehr_rare_variants = hl.import_table(DISCOVEHR_RARE_VARIANT_LIST, types={'locus':hl.tlocus(reference_genome='GRCh37'), 'alleles':hl.tarray(hl.tstr)})
#ht_discovehr_rare_variants = ht_discovehr_rare_variants.key_by(ht_discovehr_rare_variants.locus, ht_discovehr_rare_variants.alleles)

ht_initial_samples = hl.import_table(INITIAL_SAMPLES, no_header=True, key='f0')
ht_ibd_samples = hl.import_table(IBD_SAMPLES, no_header=True, key='f0')
ht_sexcheck_samples = hl.import_table(SEXCHECK_SAMPLES, no_header=True, key='f0')

mt = hl.read_matrix_table(MT)
mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))  
mt = mt.filter_cols(~hl.is_defined(ht_ibd_samples[mt.col_key]))
mt = mt.filter_cols(~hl.is_defined(ht_sexcheck_samples[mt.col_key]))
#mt = mt.filter_rows((hl.is_defined(ht_initial_rare_variants[mt.row_key]) & (~hl.is_defined(ht_gnomad_rare_variants[mt.row_key]) | ~hl.is_defined(ht_discovehr_rare_variants[mt.row_key]))) |(hl.is_defined(ht_initial_rare_variants[mt.row_key]) & hl.is_defined(ht_gnomad_rare_variants[mt.row_key]) & hl.is_defined(ht_discovehr_rare_variants[mt.row_key])))

mt = mt.filter_rows(hl.is_defined(ht_initial_rare_variants[mt.row_key]))
#mt = mt.filter_rows(hl.is_defined(ht_gnomad_rare_variants[mt.row_key]))
#mt = mt.filter_rows(hl.is_defined(ht_discovehr_rare_variants[mt.row_key]))
mt.count_rows()

mt = mt.drop('a_index', 'qual', 'rsid', 'info', 'filters', 'was_split')

sample_annotations = (hl.import_table(PHENOTYPES_TABLE).key_by('SAMPID'))
constraint_annotations = hl.read_table(ANNOTATION_TABLE)

mt = mt.annotate_cols(phenotype = sample_annotations[mt.s])
mt = mt.annotate_rows(constraint = constraint_annotations[mt.row_key])

mt = mt.annotate_cols(n_RV = hl.agg.count_where(mt.GT.is_non_ref()),
                      n_RV_SNP = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_snp(mt.alleles[0], mt.alleles[1])),
                      n_RV_indel = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_indel(mt.alleles[0], mt.alleles[1])),
                      n_coding_RV_SNP = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_snp(mt.alleles[0], mt.alleles[1]) & (mt.constraint.consequence_category != "non_coding")),
                      n_coding_RV_indel = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_indel(mt.alleles[0], mt.alleles[1]) & (mt.constraint.consequence_category != "non_coding")),
                      n_RV_PTV = hl.agg.count_where(mt.GT.is_non_ref() & (mt.constraint.consequence_category == "ptv")),
                      n_RV_damaging_missense = hl.agg.count_where(mt.GT.is_non_ref() & (mt.constraint.consequence_category == "damaging_missense")),
                      n_RV_other_missense = hl.agg.count_where(mt.GT.is_non_ref() & (mt.constraint.consequence_category == "other_missense")),
                      n_RV_synonymous = hl.agg.count_where(mt.GT.is_non_ref() & (mt.constraint.consequence_category == "synonymous")),
                      n_RV_non_coding = hl.agg.count_where(mt.GT.is_non_ref() & (mt.constraint.consequence_category == "non_coding")))

mt.cols().flatten().export(RV_FILE)
