#!/bin/bash

for i in `seq 1 22`;
do
	./plink --bfile /media/veracrypt10/Analyses/plink_prune/filterGT.chr${i} \
	--indep-pairwise 1000kb 1 0.2 --out /media/veracrypt10/Analyses/plink_prune/pruned/scz_chr${i}
done 

./plink --bfile /media/veracrypt10/Analyses/plink_prune/filterGT.chrX \
      --indep-pairwise 1000kb 1 0.2 --out /media/veracrypt10/Analyses/plink_prune/pruned/scz_chrX

cd /media/veracrypt10/Analyses/plink_prune/pruned
cat scz_chr1.prune.in scz_chr2.prune.in > scz_prune.keep.variant_list_tmp

for i in `seq 3 22`;
do
	cat scz_prune.keep.variant_list_tmp scz_chr${i}.prune.in > scz_prune.keep.variant_list
	mv scz_prune.keep.variant_list scz_prune.keep.variant_list_tmp
done

mv scz_prune.keep.variant_list_tmp scz_prune.keep.variant_list

