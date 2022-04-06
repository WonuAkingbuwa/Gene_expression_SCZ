#!/bin/bash
#SBATCH -t 00:15:00
#SBATCH -N 1

# Specify pheno names
PHENOLIST='SCZ SCZeasauto'

# Specify files
REFDATADIR=/home/akingbuw/gene_expression_project/g1000_eur
REFDATA=g1000_eur
REFDATADIR_EAS=/home/akingbuw/gene_expression_project/g1000_eas
REFDATA_EAS=g1000_eas
	# name of plink files (without extention) that will be used as reference data

# Specify output directory
OUTDIR=/home/akingbuw/gene_expression_project/Magma_output

# Specify columns names of sumstats file
SNP=SNP
PVAL=P
N=Neff
	# if no N column is present in the sumstats, the MAGMA command can be adjusted: remove "ncol=${N}" and replace it with "N=<your sample size>"


echo start of job

# copy genotype data
cp ${REFDATADIR}/${REFDATA}.bed $TMPDIR
cp ${REFDATADIR}/${REFDATA}.bim $TMPDIR
cp ${REFDATADIR}/${REFDATA}.fam $TMPDIR
cp ${REFDATADIR}/${REFDATA}.synonyms $TMPDIR

cp ${REFDATADIR_EAS}/${REFDATA_EAS}.bed $TMPDIR
cp ${REFDATADIR_EAS}/${REFDATA_EAS}.bim $TMPDIR
cp ${REFDATADIR_EAS}/${REFDATA_EAS}.fam $TMPDIR
cp ${REFDATADIR_EAS}/${REFDATA_EAS}.synonyms $TMPDIR

echo "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"
echo "GENOTYPE DATA HAS BEEN COPIED"
echo "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"

# copy sumstats files
cp /home/akingbuw/gene_expression_project/SCZeasauto_sumstats_info $TMPDIR
cp /home/akingbuw/gene_expression_project/SCZ_sumstats_info $TMPDIR

echo "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"
echo "SUMSTATS FILES HAVE BEEN COPIED"
echo "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"

# copy annotation file (MAGMA output from previous step)
for PHENO in ${PHENOLIST}; do
	(
	cp ${OUTDIR}/${PHENO}/${PHENO}.genes.annot $TMPDIR
	cd $TMPDIR
	mkdir $TMPDIR/output_${PHENO}
	) &
done
wait


cd $TMPDIR
#gunzip *.gz

# module load magma

echo "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"
echo "STARTING GENE-BASED ANALYSIS"
echo "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"

PHENO=SCZeasauto
SUMSTATS=SCZeasauto_sumstats_info
for CHR in {1..22}; do
	(
	/home/akingbuw/gene_expression_project/magma --bfile ${REFDATA_EAS} \
	--gene-annot ${PHENO}.genes.annot \
	--pval ${SUMSTATS} ncol=${N} use=${SNP},${PVAL} \
	--batch ${CHR} chr \
	--out output_${PHENO}/${PHENO}
	) &
done
wait

PHENO=SCZ
SUMSTATS=SCZ_sumstats_info
for CHR in {1..22}; do
	(
	/home/akingbuw/gene_expression_project/magma --bfile ${REFDATA} \
	--gene-annot ${PHENO}.genes.annot \
	--pval ${SUMSTATS} N=105318 use=${SNP},${PVAL} \
	--batch ${CHR} chr \
	--out output_${PHENO}/${PHENO}
	) &
done
wait

#echo "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"
#echo "WHAT'S IN TMPDIR/output?"
#ls $TMPDIR/output
#echo "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"



for PHENO in ${PHENOLIST}; do
	(
	cp $TMPDIR/output_${PHENO}/* $OUTDIR/${PHENO}
	) &
done
wait

echo end of job
