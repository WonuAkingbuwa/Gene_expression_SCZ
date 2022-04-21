library(ggplot2)
library(ggsci)
library(dplyr)
library(data.table)

source("Scripts/r_options_Wonu.R")

df_after <- fread(SAMPLE_AFTER_QC_FILE, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE) %>%
  mutate(phase='After Variant QC')

#read in pheno file 
df_pheno <- fread('data/scz_phenoinfo.txt', stringsAsFactors=FALSE, sep='\t', header=TRUE, data.table=FALSE) %>% rename(s=SAMPID)

if(any(duplicated(df_pheno$s))) {
  df_pheno <- df_pheno[-which(duplicated(df_pheno$s)),]
}

df_after_merged <- merge(df_after, df_pheno, by='s')

df_keep <- df_after['s']
print(paste0("Started with: ", nrow(df_keep), " samples")) #12070

n_MADs <- 4

# r_ti_tv
df_keep_ti_tv <- group_by(df_after_merged, SAMP_SOURCE) %>%
  summarise(median=median(sample_qc.r_ti_tv), mad=mad(sample_qc.r_ti_tv)) %>%
  inner_join(df_after_merged, by='SAMP_SOURCE') %>%
  filter((sample_qc.r_ti_tv >= median - n_MADs*mad & sample_qc.r_ti_tv <= median + n_MADs*mad) | is.na(mad))

df_keep <- df_keep_ti_tv %>% inner_join(df_keep, by='s')
print(paste0("Remove Ti/Tv outliers: ", nrow(df_keep), " samples remain")) 

# r_het_hom_var
df_keep_het_hom_var <- group_by(df_after_merged, SAMP_SOURCE) %>%
  summarise(median=median(sample_qc.r_het_hom_var), mad=mad(sample_qc.r_het_hom_var)) %>%
  inner_join(df_after_merged, by='SAMP_SOURCE') %>%
  filter((sample_qc.r_het_hom_var >= median - n_MADs*mad & sample_qc.r_het_hom_var <= median + n_MADs*mad) | is.na(mad))

df_keep <- df_keep_het_hom_var %>% inner_join(df_keep, by='s')
print(paste0("Remove Het/HomVar outliers: ", nrow(df_keep), " samples remain")) 

# r_insertion_deletion
df_keep_insertion_deletion <- group_by(df_after_merged, SAMP_SOURCE) %>%
  summarise(median=median(sample_qc.r_insertion_deletion), mad=mad(sample_qc.r_insertion_deletion)) %>%
  inner_join(df_after_merged, by='SAMP_SOURCE') %>%
  filter((sample_qc.r_insertion_deletion >= median - n_MADs*mad & sample_qc.r_insertion_deletion <= median + n_MADs*mad) | is.na(mad))

df_keep <- df_keep_insertion_deletion %>% inner_join(df_keep, by='s')
print(paste0("Remove Ins/Del outliers: ", nrow(df_keep), " samples remain")) 

# n_snps
df_keep_n_snps <- group_by(df_after_merged, SAMP_SOURCE) %>%
  summarise(median=median(sample_qc.n_snp), mad=mad(sample_qc.n_snp)) %>%
  inner_join(df_after_merged, by='SAMP_SOURCE') %>%
  filter((sample_qc.n_snp >= median - n_MADs*mad & sample_qc.n_snp <= median + n_MADs*mad) | is.na(mad))

df_keep <- df_keep_n_snps %>% inner_join(df_keep, by='s')
print(paste0("Remove n_snp outliers: ", nrow(df_keep), " samples remain")) 

#n_insertions
# n_insertion
df_keep_insertion <- group_by(df_after_merged, SAMP_SOURCE) %>%
  summarise(median=median(sample_qc.n_insertion), mad=mad(sample_qc.n_insertion)) %>%
  inner_join(df_after_merged, by='SAMP_SOURCE') %>%
  filter((sample_qc.n_insertion >= median - n_MADs*mad & sample_qc.n_insertion <= median + n_MADs*mad) | is.na(mad))

df_keep <- df_keep_insertion %>% inner_join(df_keep, by='s')
print(paste0("Remove n_insertion outliers: ", nrow(df_keep), " samples remain")) 

#n_deletions
# n_deletion
df_keep_deletion <- group_by(df_after_merged, SAMP_SOURCE) %>%
  summarise(median=median(sample_qc.n_deletion), mad=mad(sample_qc.n_deletion)) %>%
  inner_join(df_after_merged, by='SAMP_SOURCE') %>%
  filter((sample_qc.n_deletion >= median - n_MADs*mad & sample_qc.n_deletion <= median + n_MADs*mad) | is.na(mad))

df_keep <- df_keep_deletion %>% inner_join(df_keep, by='s')
print(paste0("Remove n_deletion outliers: ", nrow(df_keep), " samples remain")) 

df_final_sample_summary <- data.table(Filter = c("Samples after population filters",
                          paste0("Within batch Ti/Tv ratio outside ", n_MADs, " median deviations"),
                          paste0("Within batch Het/HomVar ratio outside ", n_MADs, " median deviations"),
                          paste0("Within batch Insertion/Deletion ratio outside ", n_MADs, " median deviations"),
                          paste0("Within location n snps outside ", n_MADs, " median deviations"),
                          paste0("Within location n insertions outside ", n_MADs, " median deviations"),
                          paste0("Within location n deletions outside ", n_MADs, " median deviations"),
                          "Samples after final sample filters"),
                     "Samples" = c(nrow(df_after_merged),
                                nrow(df_after_merged) - nrow(df_keep_ti_tv),
                                nrow(df_after_merged) - nrow(df_keep_het_hom_var),
                                nrow(df_after_merged) - nrow(df_keep_insertion_deletion),
                                nrow(df_after_merged) - nrow(df_keep_n_snps),
                                nrow(df_after_merged) - nrow(df_keep_insertion),
                                nrow(df_after_merged) - nrow(df_keep_deletion),
                                nrow(df_keep)),
                     'Cases' = c(nrow(df_after_merged %>% filter(ANALYSIS_CAT=="Case")),
                      nrow(df_after_merged %>% filter(ANALYSIS_CAT=="Case")) - nrow(df_keep_ti_tv %>% filter(ANALYSIS_CAT == "Case")),
                      nrow(df_after_merged %>% filter(ANALYSIS_CAT=="Case")) - nrow(df_keep_het_hom_var %>% filter(ANALYSIS_CAT == "Case")),
                      nrow(df_after_merged %>% filter(ANALYSIS_CAT=="Case")) - nrow(df_keep_insertion_deletion %>% filter(ANALYSIS_CAT == "Case")),
                      nrow(df_after_merged %>% filter(ANALYSIS_CAT=="Case")) - nrow(df_keep_n_snps %>% filter(ANALYSIS_CAT == "Case")),
                      nrow(df_after_merged %>% filter(ANALYSIS_CAT=="Case")) - nrow(df_keep_insertion %>% filter(ANALYSIS_CAT == "Case")),
                      nrow(df_after_merged %>% filter(ANALYSIS_CAT=="Case")) - nrow(df_keep_deletion %>% filter(ANALYSIS_CAT == "Case")),
                      nrow(df_keep %>% filter(ANALYSIS_CAT.x=="Case"))),
                     "Controls" = c(nrow(df_after_merged %>% filter(ANALYSIS_CAT=="Control")),
                      nrow(df_after_merged %>% filter(ANALYSIS_CAT=="Control")) - nrow(df_keep_ti_tv %>% filter(ANALYSIS_CAT == "Control")),
                      nrow(df_after_merged %>% filter(ANALYSIS_CAT=="Control")) - nrow(df_keep_het_hom_var %>% filter(ANALYSIS_CAT == "Control")),
                      nrow(df_after_merged %>% filter(ANALYSIS_CAT=="Control")) - nrow(df_keep_insertion_deletion %>% filter(ANALYSIS_CAT == "Control")),
                      nrow(df_after_merged %>% filter(ANALYSIS_CAT=="Control")) - nrow(df_keep_n_snps %>% filter(ANALYSIS_CAT == "Control")),
                      nrow(df_after_merged %>% filter(ANALYSIS_CAT=="Control")) - nrow(df_keep_insertion %>% filter(ANALYSIS_CAT == "Control")),
                      nrow(df_after_merged %>% filter(ANALYSIS_CAT=="Control")) - nrow(df_keep_deletion %>% filter(ANALYSIS_CAT == "Control")),
                      nrow(df_keep %>% filter(ANALYSIS_CAT.x=="Control"))))

fwrite(df_final_sample_summary, file='Analyses/15_sample_count.tsv', quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')

# write out
df_keep <- df_keep %>% select(s)
print(dim(df_keep))

fwrite(df_keep, file=FINAL_SAMPLE_LIST, quote=FALSE, row.names=FALSE, col.names=FALSE)
