library(data.table)
library(ggplot2)
library(reshape2)

#rare PTVs####
#import geneset regression estimates (total + exclusions) + add data identifiers to the header names
ptv_tot_estimates <- fread("scz_rv_ptv_geneset_estimates_noURV_ewptvcov.tsv", sep='\t', header=TRUE, data.table=FALSE)
ptv_tot_estimates$geneset <- gsub("n_RV_PTV__RVs_|_intervals_bed_tsv", "", ptv_tot_estimates$geneset)
ptv_tot_estimates <- ptv_tot_estimates[,c(2:9)]
colnames(ptv_tot_estimates)[2:8] <- paste("tot", colnames(ptv_tot_estimates[,c(2:8)]), sep = "_")

ptv_FIN_estimates <- fread("scz_rv_ptv_geneset_estimates_FIN_exclusions.tsv", sep='\t', header=TRUE, data.table=FALSE)
ptv_FIN_estimates$geneset <- gsub("n_RV_PTV__RVs_|_intervals_bed_tsv", "", ptv_FIN_estimates$geneset)
ptv_FIN_estimates <- ptv_FIN_estimates[,c(2:9)]
colnames(ptv_FIN_estimates)[2:8] <- paste("fin", colnames(ptv_FIN_estimates[,c(2:8)]), sep = "_")

ptv_NorSwe_estimates <- fread("scz_rv_ptv_geneset_estimates_NorSwe_exclusions.tsv", sep='\t', header=TRUE, data.table=FALSE)
ptv_NorSwe_estimates$geneset <- gsub("n_RV_PTV__RVs_|_intervals_bed_tsv", "", ptv_NorSwe_estimates$geneset)
ptv_NorSwe_estimates <- ptv_NorSwe_estimates[,c(2:9)]
colnames(ptv_NorSwe_estimates)[2:8] <- paste("nor", colnames(ptv_NorSwe_estimates[,c(2:8)]), sep = "_")

ptv_othEUR_estimates <- fread("scz_rv_ptv_geneset_estimates_othEUR_exclusions.tsv", sep='\t', header=TRUE, data.table=FALSE)
ptv_othEUR_estimates$geneset <- gsub("n_RV_PTV__RVs_|_intervals_bed_tsv", "", ptv_othEUR_estimates$geneset)
ptv_othEUR_estimates <- ptv_othEUR_estimates[,c(2:9)]
colnames(ptv_othEUR_estimates)[2:8] <- paste("oth", colnames(ptv_othEUR_estimates[,c(2:8)]), sep = "_")

ptv_allEURout_estimates <- fread("scz_rv_ptv_geneset_estimates_allEURout_exclusions.tsv", sep='\t', header=TRUE, data.table=FALSE)
ptv_allEURout_estimates$geneset <- gsub("n_RV_PTV__RVs_|_intervals_bed_tsv", "", ptv_allEURout_estimates$geneset)
ptv_allEURout_estimates <- ptv_allEURout_estimates[,c(2:9)]
colnames(ptv_allEURout_estimates)[2:8] <- paste("all", colnames(ptv_allEURout_estimates[,c(2:8)]), sep = "_")

#merge data
all_estimates <- merge(ptv_tot_estimates, ptv_FIN_estimates, by = "geneset")
all_estimates <- merge(all_estimates, ptv_NorSwe_estimates, by = "geneset")
all_estimates <- merge(all_estimates, ptv_othEUR_estimates, by = "geneset")
all_estimates <- merge(all_estimates, ptv_allEURout_estimates, by = "geneset")
all_estimates <- all_estimates[,c(1,2,9,16,23,30)]

#plot
tiff("pca_exclusions_plot_rv_ptv.tiff", width = 3000, height = 2000, res=600)
  ggplot(all_estimates, aes(tot_beta)) +
  geom_point(aes(y=fin_beta, colour="Finnish"), size=0.5) +
  geom_point(aes(y=nor_beta, colour="Northern Swedish"), size=0.5) +
  geom_point(aes(y=oth_beta, colour="Other European"), size=0.5) +
  geom_point(aes(y=all_beta, colour="All outliers"), size=0.5) +
  
  ylab("Beta (Samples excluding outlier populations)") +
  xlab("Beta (total European sample)") +
  theme(axis.title.x = element_text(size=6)) +
  theme(axis.title.y = element_text(size=6)) +
  theme(axis.text.x = element_text(size=5)) +
  theme(axis.text.y = element_text(size=5)) +
    
  labs(colour="Populations excluded") +
  theme(legend.text = element_text(size = 4)) +
  theme(legend.title = element_text(size = 5))
    
  dev.off()
  
