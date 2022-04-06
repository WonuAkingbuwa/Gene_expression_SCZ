#!/bin/bash
#SBATCH -t 00:05:00
#SBATCH -N 1

DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=complete_genesets_transpose.gmt
PHENOLIST='SCZ SCZeasauto'

cp ${DIR}/${SET} $TMPDIR

for PHENO in ${PHENOLIST}; do
	(
	cp ${DIR}/${PHENO}/${PHENO}.genes.raw $TMPDIR
	cd $TMPDIR
	mkdir $TMPDIR/output_${PHENO}
	) &
done
wait

cd $TMPDIR
# module load magma


# Run gene-set analysis

for PHENO in ${PHENOLIST}; do
	(
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait
