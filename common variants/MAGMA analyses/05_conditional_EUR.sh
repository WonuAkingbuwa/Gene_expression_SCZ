#!/bin/bash
#SBATCH -t 02:00:00
#SBATCH -N 1

DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_pipyca1_both.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=pyramidal_CA1,PI_genes --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_pipyca1_pi.gmt.txt
PHENOLIST='SCZ'

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
SET=genesets_transpose_conditional_pipyca1_pyca1.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=pyramidal_CA1 --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_pipyss_both.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=pyramidal_SS,PI_genes --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_pipyss_pi.gmt.txt
PHENOLIST='SCZ'

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
SET=genesets_transpose_conditional_pipyss_pyss.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=pyramidal_SS --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_piserneu_both.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=Serotonergic_Neuron,PI_genes --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_piserneu_pi.gmt.txt
PHENOLIST='SCZ'

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
SET=genesets_transpose_conditional_piserneu_serneu.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=Serotonergic_Neuron --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_pistrint_both.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=Striatal_Interneuron,PI_genes --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_pistrint_pi.gmt.txt
PHENOLIST='SCZ'

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
SET=genesets_transpose_conditional_pistrint_strint.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=Striatal_Interneuron --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_pida_both.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=Dopaminergic_Adult,PI_genes --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_pida_da.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=Dopaminergic_Adult --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_pida_pi.gmt.txt
PHENOLIST='SCZ'

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
SET=genesets_transpose_conditional_piexca3_exca3.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=exCA3 --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_piexca3_pi.gmt.txt
PHENOLIST='SCZ'

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
PHENOLIST='SCZ'

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
PHENOLIST='SCZ'

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
PHENOLIST='SCZ'

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
SET=genesets_transpose_conditional_piexpfc1_both.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=exPFC1,PI_genes --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_piexpfc1_expfc1.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=exPFC1 --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_piexpfc1_pi.gmt.txt
PHENOLIST='SCZ'

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
SET=genesets_transpose_conditional_piexpfc2_both.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=exPFC2,PI_genes --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_piexpfc2_expfc2.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=exPFC2 --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_piexpfc2_pi.gmt.txt
PHENOLIST='SCZ'

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
SET=genesets_transpose_conditional_pigaba1_both.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=GABA1,PI_genes --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_pigaba1_gaba1.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=GABA1 --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_pigaba1_pi.gmt.txt
PHENOLIST='SCZ'

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
SET=genesets_transpose_conditional_pigaba2_both.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=GABA2,PI_genes --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_pigaba2_gaba2.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=GABA2 --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_pigaba2_pi.gmt.txt
PHENOLIST='SCZ'

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
SET=genesets_transpose_conditional_piint_both.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=interneurons,PI_genes --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_piint_int.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=interneurons --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_piint_pi.gmt.txt
PHENOLIST='SCZ'

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
SET=genesets_transpose_conditional_pimg_both.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=MG,PI_genes --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_pimg_mg.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=MG --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_pimg_pi.gmt.txt
PHENOLIST='SCZ'

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
SET=genesets_transpose_conditional_pimsn_both.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=Medium_Spiny_Neuron,PI_genes --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_pimsn_msn.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=Medium_Spiny_Neuron --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_pimsn_pi.gmt.txt
PHENOLIST='SCZ'

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
PHENOLIST='SCZ'

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
PHENOLIST='SCZ'

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
PHENOLIST='SCZ'

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
SET=genesets_transpose_conditional_piodc1_both.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=ODC1,PI_genes --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_piodc1_odc1.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=ODC1 --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_piodc1_pi.gmt.txt
PHENOLIST='SCZ'

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
SET=genesets_transpose_conditional_piopc_both.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=OPC,PI_genes --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_piopc_opc.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=OPC --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_piopc_pi.gmt.txt
PHENOLIST='SCZ'

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
PHENOLIST='SCZ'

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


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_pidn_dn.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=Dopaminergic_Neuroblast --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_pidn_pi.gmt.txt
PHENOLIST='SCZ'

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
SET=genesets_transpose_conditional_piedn_both.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=Embryonic_Dopaminergic_Neuron,PI_genes --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_piedn_edn.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=Embryonic_Dopaminergic_Neuron --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_piedn_pi.gmt.txt
PHENOLIST='SCZ'

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
SET=genesets_transpose_conditional_piexca1_both.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=exCA1,PI_genes --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_piexca1_exca1.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=exCA1 --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait


DIR=/home/akingbuw/gene_expression_project/Magma_output
SET=genesets_transpose_conditional_piexca1_pi.gmt.txt
PHENOLIST='SCZ'

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
SET=genesets_transpose_conditional_piexca3_both.gmt.txt
PHENOLIST='SCZ'

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
	/home/akingbuw/gene_expression_project/magma --gene-results ${PHENO}.genes.raw --set-annot ${DIR}/${SET} --model condition-hide=exCA3,PI_genes --out output_${PHENO}/${PHENO}_${SET}
	cp $TMPDIR/output_${PHENO}/* ${DIR}/${PHENO}
	) &
done
wait
