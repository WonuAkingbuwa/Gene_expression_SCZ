#!/bin/bash
#SBATCH -t 02:00:00
#SBATCH -N 1

DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_pidn_dn.gmt.txt
PHENOLIST='SCZeasauto'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=Dopaminergic_Neuroblast --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_pidn_pi.gmt.txt
PHENOLIST='SCZeasauto'

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



for PHENO in ${PHENOLIST}; do
	(
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=PI_genes --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_pineur_both.gmt.txt
PHENOLIST='SCZeasauto'

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



for PHENO in ${PHENOLIST}; do
	(
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=Neuroblasts,PI_genes --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_pineur_neur.gmt.txt
PHENOLIST='SCZeasauto'

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



for PHENO in ${PHENOLIST}; do
	(
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=Neuroblasts --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_pineur_pi.gmt.txt
PHENOLIST='SCZeasauto'


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



for PHENO in ${PHENOLIST}; do
	(
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=PI_genes --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_piexdg_both.gmt.txt
PHENOLIST='SCZeasauto'

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



for PHENO in ${PHENOLIST}; do
	(
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=exDG,PI_genes --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_piexdg_exdg.gmt.txt
PHENOLIST='SCZeasauto'

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



for PHENO in ${PHENOLIST}; do
	(
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=exDG --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_piexdg_pi.gmt.txt
PHENOLIST='SCZeasauto'

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



for PHENO in ${PHENOLIST}; do
	(
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=PI_genes --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_pidn_both.gmt.txt
PHENOLIST='SCZeasauto'

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



for PHENO in ${PHENOLIST}; do
	(
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=Dopaminergic_Neuroblast,PI_genes --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait
