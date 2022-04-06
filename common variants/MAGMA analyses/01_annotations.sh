#!/bin/bash
#SBATCH -t 00:05:00
#SBATCH -N 1

# Specify pheno names
PHENOLIST='SCZ SCZeasauto'

# Specify geneloc files
GENELOCDIR=/home/akingbuw/gene_expression_project/gene_loc_files
GENELOC=geneloc

# Specify output directory
OUTDIR=/home/akingbuw/gene_expression_project/Magma_output



echo start of job

################################
# Copy files and prepare for analysis

cp ${GENELOCDIR}/${GENELOC} $TMPDIR
cd $TMPDIR
#gunzip ${GENELOC}.gz

for PHENO in ${PHENOLIST}; do
	(
	cp /home/akingbuw/gene_expression_project/${PHENO}.snp.loc $TMPDIR
	cd $TMPDIR
	#gunzip ${PHENO}.snp.loc.gz
	mkdir $TMPDIR/output_${PHENO}

	) &
done
wait

cd $TMPDIR


################################
# Annotate SNPs to genes
# Copy files back to home dir

for PHENO in ${PHENOLIST}; do
	(
	/home/akingbuw/gene_expression_project/magma \
	--annotate \
	--snp-loc ${PHENO}.snp.loc \
	--gene-loc ${GENELOC} \
	--out output_${PHENO}/${PHENO}

	cp $TMPDIR/output_${PHENO}/* $OUTDIR/${PHENO}
	) &
done
wait


echo end of job
