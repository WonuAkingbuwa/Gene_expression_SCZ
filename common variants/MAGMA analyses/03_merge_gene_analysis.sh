#run this in the folder with all the previous output
cd Magma_output

module load magma

PHENO=SCZeasauto
/home/akingbuw/gene_expression_project/magma --merge ${PHENO}/${PHENO} --out ${PHENO}/${PHENO}

PHENO=SCZ
/home/akingbuw/gene_expression_project/magma --merge ${PHENO}/${PHENO} --out ${PHENO}/${PHENO}
