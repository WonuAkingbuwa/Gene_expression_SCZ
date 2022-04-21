#!/bin/bash

for i in `seq 1 22`;
do
	./plink --bfile /media/veracrypt10/Analyses/plink_prune/final_qc.chr${i} \
	--indep-pairwise 1000kb 1 0.2 --out /media/veracrypt10/Analyses/plink_prune/pruned/final_qc.chr${i}
done 

./plink --bfile /media/veracrypt10/Analyses/plink_prune/final_qc.chrX \
      --indep-pairwise 1000kb 1 0.2 --out /media/veracrypt10/Analyses/plink_prune/pruned/final_qc.chrX

cd /media/veracrypt10/Analyses/plink_prune/pruned
cat final_qc.chr1.prune.in final_qc.chr2.prune.in > prune.final_qc.keep.variant_list_tmp

for i in `seq 3 22`;
do
	cat prune.final_qc.keep.variant_list_tmp final_qc.chr${i}.prune.in > prune.final_qc.keep.variant_list
	mv prune.final_qc.keep.variant_list prune.final_qc.keep.variant_list_tmp

done

mv prune.final_qc.keep.variant_list_tmp prune.final_qc.keep.variant_list
cp prune.final_qc.keep.variant_list /media/veracrypt10/Analyses/Analyses

