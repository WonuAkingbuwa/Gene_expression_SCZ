####Import hail####
import hail as hl
hl.init(min_block_size=128)

#import other python libraries that will be relevant later
from hail.plot import show
from pprint import pprint
hl.plot.output_notebook() 

#import vcf files and write as mt for working in hail

#swedish schizophrenia data
hl.import_vcf('/media/veracrypt10/dbGAP/SWE_SCZ/PhenoGenotypeFiles/RootStudyConsentSet_phs000473.SchizophreniaSwedish_Sklar.v2.p2.c1.GRU/GenotypeFiles/phg000773.v1.SchizophreniaSwedish_Sklar_v2.genotype-calls-vcf.WES_markerset_grc37.c1.GRU/sub/sub20160809/swe.vcf').write('/media/veracrypt10/Analyses/Analyses/scz_main.mt', overwrite=True)

#1KG data for PCA analyses down the line
hl.import_vcf('/media/veracrypt10/Analyses/data/*.vcf.bgz').write('Analyses/1KG_autosomes_GRCh37.mt', overwrite=True)
