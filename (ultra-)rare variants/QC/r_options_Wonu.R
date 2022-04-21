# 03_initial_sample_qc_filter.r
QC_FILE <- "Analyses/initial_sample_qc.tsv"
SAMPLE_LIST_INITIAL_QC <- 'Analyses/initial_qc.keep.sample_list'

# 03_initial_sample_qc_plot.r 
PLOTS <- 'Output/'
# Define some thresholds 
T_sample_callRate <- 0
T_dpMean <- 0
T_gqMean <- 0

# 05_impute_sex_plot.r
IMPUTESEX_FILE <- "Analyses/imputesex.tsv"
SEXCHECK_LIST <- 'Analyses/sexcheck.remove.sample_list'
Y_NCALLED_FILE <- "Analyses/ycalled.tsv"
T_impute_sex <- 0.6

# 06_ibd_plot.r
IBD_FILE <- "Analyses/ibd.tsv"
IBD_THRESHOLD <- 0.2

# 06_ibd_filtered.r
SAMPLE_LIST_IBD <- 'Analyses/ibd.remove.sample_list'

# 08_ultra_rare_counts_plot.r 
URV_FILE <- 'annotations/URVs.tsv'
  
# 09_10_pca_plot.r
PCA_SCORES <- 'Analyses/pca_scores.tsv'
PCA_1KG_SCORES <- 'Analyses/pca_scores_1kg.tsv'
EUROPEAN_SAMPLES_STRICT <- 'Analyses/european.strict.sample_list'
EUROPEAN_SAMPLES_LOOSE <- 'Analyses/european.loose.sample_list'
EUROPEAN_SAMPLES_EXCLUDING_URV_OUTLIERS <- 'Analyses/european.no_URV_outliers.sample_list'

T_nURVSNP <- 300
T_nURVIndel <- 25
  
T_European_RF <- 0.95

PCA_EUR_SCORES <- "Analyses/pca_scores.strict_european.tsv"

# 13_final_variant_qc_plot.r
VARIANT_QC_FILE <- 'Analyses/final_qc.variants.tsv.gz'
T_variant_call_rate  <- 0.80
T_absdiff <- 0.02
T_pHWE <- 1e-6

# 13_final_variant_qc_filter.r
VARIANT_LIST <- 'Analyses/final_qc.keep.variant_list'

# 14_final_sample_qc_plot.r
SAMPLE_BEFORE_QC_FILE <- 'Analyses/final_qc.before.samples.tsv'
SAMPLE_AFTER_QC_FILE <- 'Analyses/final_qc.after.samples.tsv'

# 14_final_sample_qc_filter.r
FINAL_SAMPLE_LIST <- 'Analyses/final_qc.keep.sample_list'
FINAL_SAMPLE_LIST_REMOVE_SINGLETON_OUTLIERS <- 'Analyses/final_qc_remove_singleton_outliers.keep.sample_list'