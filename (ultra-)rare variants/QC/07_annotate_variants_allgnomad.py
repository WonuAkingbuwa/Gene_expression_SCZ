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

from gnomad.utils.vep import process_consequences

####import/download gnomad, discoverHER, and mpc files for variant annotation

hl.import_vcf('/media/veracrypt10/Analyses/annotations/GHS_Freeze_50.L3DP10.pVCF.frq.vcf', reference_genome='GRCh37').write('annotations/discovEHR.mt', overwrite=True)
DISCOVEHR_SITES_MT = 'annotations/discovEHR.mt'
discovEHR_mt = hl.read_matrix_table(DISCOVEHR_SITES_MT)
discovEHR_mt.describe()
discovEHR_mt.count_rows() 
discovEHR_mt = discovEHR_mt.filter_rows(hl.len(discovEHR_mt.filters) == 0)
discovEHR_mt.count_rows() 
discovEHR_ht = discovEHR_mt.rows()
discovEHR_ht.describe()
discovEHR_ht.write("annotations/discovEHR_variants.ht", overwrite=True)

#import MPC (missense badness) scores and write out to ht for future use - currently mapped to GRCh37, no need for liftover
#MPC_SCORE = 'annotations/fordist_constraint_official_mpc_values_v2.txt.bgz'
#mpc_ht = hl.import_table(MPC_SCORE, impute=True)
#mpc_ht.describe()
#mpc_ht = mpc_ht.annotate(locus = hl.locus(contig = mpc_ht.chrom, pos = mpc_ht.pos),
#    alleles = [mpc_ht.ref, mpc_ht.alt])
#mpc_ht.describe()
#mpc_ht = mpc_ht.key_by(mpc_ht.locus, mpc_ht.alleles).select('MPC')
#mpc_ht.describe()
#mpc_ht.write('annotations/fordist_constraint_official_mpc_values_v2_GRCh37.ht', overwrite=True)

####annotate variants with VEP#### - consider filtering to rare variants here 
#import all entries filtered file 
MT = 'Analyses/scz_allentries.mt'
ht = hl.read_matrix_table(MT).rows()

#run VEP annotation command and write out VEP ht
ht_vep = hl.vep(ht, "Scripts/vep_configuration_wonu_020720.json")
ht_vep.write("annotations/vep_annotate_wonu.ht", overwrite=True)
ht_vep.write("annotations/vep_annotate_wonu_norarefilter.ht", overwrite=True)

#annotation file to be written out 
ANNOTATION_TABLE = 'annotations/gene.ht'

#read in all annotation files
GNOMAD_SITES_HT = 'annotations/gnomad.exomes.r2.1.1.sites.ht'
gnomad_ht = hl.read_table(GNOMAD_SITES_HT)
gnomad_ht.describe()
gnomad_ht.freq.AC[0].show() 

#filter AC < 2
gnomad_ht_AC = gnomad_ht.filter((gnomad_ht.freq.AC[0] == 0) | (gnomad_ht.freq.AC[0] == 1))
gnomad_ht_AC.freq.AC[0].show()

#GNOMAD_SITES_MT = 'annotations/gnomad.exomes.r2.1.1.sites.mt'
#gnomad_mt = hl.read_matrix_table(GNOMAD_SITES_MT)
#gnomad_mt.describe()
#gnomad_mt.info.AC.show() 

DISCOVHER_VARIANTS = 'annotations/discovEHR_variants.ht'

MPC_SCORE = 'annotations/fordist_constraint_official_mpc_values_v2_GRCh37.ht'
mpc_ht = hl.read_table(MPC_SCORE)
mpc_ht = mpc_ht.select(mpc_ht.MPC)

mt = hl.read_matrix_table(MT)

# Annotate with the vep information
vep_ht = hl.read_table("annotations/vep_annotate_wonu_norarefilter.ht")
vep_ht.describe()
vep_ht = vep_ht.select('vep')
mt = mt.annotate_rows(vep_ann = vep_ht[mt.row_key])
mt = mt.annotate_rows(vep = mt.vep_ann.vep)
mt = mt.drop(mt.vep_ann)

#annotate with MPC info
mt = mt.annotate_rows(mpc = mpc_ht[mt.row_key])

#specify annotations for the different variant classes
ptv = hl.set(["transcript_ablation", "splice_acceptor_variant",
              "splice_donor_variant", "stop_gained", "frameshift_variant"])

missense = hl.set(["stop_lost", "start_lost", "transcript_amplification",
                   "inframe_insertion", "inframe_deletion", "missense_variant",
                   "protein_altering_variant", "splice_region_variant"])

synonymous = hl.set(["incomplete_terminal_codon_variant", "stop_retained_variant", "synonymous_variant"])

non_coding = hl.set(["coding_sequence_variant", "mature_miRNA_variant", "5_prime_UTR_variant",
              "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
              "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant",
              "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
              "regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation",
              "regulatory_region_variant", "feature_truncation", "intergenic_variant"])

#add annotation for variant and gene/consequence.
mt = process_consequences(mt) #this is from the gnomad utils function

mt = mt.annotate_rows(consequence_category = 
    hl.case().when(ptv.contains(mt.vep.worst_csq_for_variant_canonical.most_severe_consequence), "ptv")
             .when(missense.contains(mt.vep.worst_csq_for_variant_canonical.most_severe_consequence) & 
                   (~hl.is_defined(mt.vep.worst_csq_for_variant_canonical.polyphen_prediction) | 
                    ~hl.is_defined(mt.vep.worst_csq_for_variant_canonical.sift_prediction) ), "other_missense")
             .when(missense.contains(mt.vep.worst_csq_for_variant_canonical.most_severe_consequence) & 
                   (mt.vep.worst_csq_for_variant_canonical.sift_pred == "D") & 
                   (mt.vep.worst_csq_for_variant_canonical.polyphen2_hdiv_pred == "D") & 
                   (mt.vep.worst_csq_for_variant_canonical.polyphen2_hvar_pred == "D") & 
                   (mt.vep.worst_csq_for_variant_canonical.lrt_pred == "D") & 
                   ((mt.vep.worst_csq_for_variant_canonical.mutationtaster_pred == "A") | (mt.vep.worst_csq_for_variant_canonical.mutationtaster_pred == "D")) & 
                   ((mt.vep.worst_csq_for_variant_canonical.mutationassessor_pred == "H") | (mt.vep.worst_csq_for_variant_canonical.mutationassessor_pred == "M")) & 
                   (mt.vep.worst_csq_for_variant_canonical.provean_pred == "D"), "damaging_missense")
             .when(missense.contains(mt.vep.worst_csq_for_variant_canonical.most_severe_consequence), "other_missense")
             .when(synonymous.contains(mt.vep.worst_csq_for_variant_canonical.most_severe_consequence), "synonymous")
             .when(non_coding.contains(mt.vep.worst_csq_for_variant_canonical.most_severe_consequence), "non_coding")
             .default("NA")
    )

# Now that we have the gene information for each variant, we can evaluate the constraint.
mt = mt.annotate_rows(inGnomAD = hl.is_defined(gnomad_ht_AC[mt.row_key]))

discovehr_variants_ht = hl.read_table(DISCOVHER_VARIANTS)
mt = mt.annotate_rows(indiscovEHR = hl.is_defined(discovehr_variants_ht[mt.row_key]))

#mt = mt.checkpoint('annotations/annotations_checkpoint.mt')
#mt = hl.read_matrix_table('annotations/annotations_checkpoint.mt')

mt_rows = mt.rows().repartition(64)
mt_rows.write(ANNOTATION_TABLE, overwrite=True) #if this does not work try checkpoint instead of write. if this fails too, try checkpoint command after annotating with gnomad and discovehr 
















