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

###import raw hail mt and prep for quality control###
mt = hl.read_matrix_table('Analyses/scz_main.mt') 

# Count before splitting multi-allelics.
n = mt.count()

pprint('n samples:')
print(n[1]) 
pprint('n variants:')
print(n[0])

#filter out variants with more than 6 alleles###
mt = mt.filter_rows(mt.alleles.length() <= 6)
n = mt.count_rows()
pprint('')
pprint('n variants not more than 6 alleles:')
print(n)

#split multi allelic variants### - adding "permit_shuffle=True" fixes downstream duplicate loci error, at the cost of additional compute time
mt = hl.split_multi_hts(mt, permit_shuffle=True)

mt = mt.filter_entries(
    hl.is_defined(mt.GT) &
    (
        (mt.GT.is_hom_ref() & 
            (
               # ((mt.AD[0] / mt.DP) < 0.8) | 
                (mt.GQ < 20) |
                (mt.DP < 10)
        	)
        ) |
        (mt.GT.is_het() & 
        	( 
                (((mt.AD[0] + mt.AD[1]) / mt.DP) < 0.8) | 
                ((mt.AD[1] / mt.DP) < 0.2) | 
                (mt.PL[0] < 20) |
                (mt.DP < 10)
        	)
        ) |
        (mt.GT.is_hom_var() & 
        	(
                ((mt.AD[1] / mt.DP) < 0.8) |
                (mt.PL[0] < 20) |
                (mt.DP < 10)
        	)
        )
    ),
    keep = False
)

#write out filtered full entries and hardcalls files
mt.write('Analyses/scz_allentries.mt', overwrite = True)
mt = mt.checkpoint('Analyses/scz_allentries.mt', overwrite = True)
mt = hl.read_matrix_table('Analyses/scz_allentries.mt')

mt.select_entries(mt.GT).repartition(100).write('Analyses/scz_hardcalls.mt', overwrite = True)

mt = hl.read_matrix_table('Analyses/scz_hardcalls.mt')
n = mt.count()

pprint('n samples:')
print(n[1]) 
pprint('n variants:')
print(n[0]) 
