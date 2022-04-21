library(dplyr)
library(data.table)
source("r_options_Wonu.R")

# Run the plotting again to ensure that the thresholds are as in the plots.
source("sample_qc_plots.R")

QC_FILE <- "Analyses/initial_sample_qc.tsv"
df <- fread(QC_FILE, sep='\t', header=TRUE, data.table=FALSE)
names(df) <- gsub("sample_qc\\.", "", names(df))
names(df) <- gsub("phenotype\\.", "", names(df))

if(any(duplicated(df$s))) {
  cat('remove duplicated names...\n')
  df <- df[-which(duplicated(df$s)),]
}

df_merged <- df

df_initial_summary_count <- data.table(
  Filter = c(
    "Initial samples in vcf", 
    "Unable to phenotype and sequence information",
    "Unknown phenotype"),
  Samples = c(nrow(df),
              nrow(df) - nrow(df_merged),
              nrow(filter(df_merged, PRIMARY_DISEASE=='Unknown'))),
  "Schizophrenia cases" = c(nrow(df_merged %>% filter(PRIMARY_DISEASE == "Schizophrenia")),
                      NA,
                      NA),
  "Bipolar cases" = c(nrow(df_merged %>% filter(PRIMARY_DISEASE == "Bipolar_Disorder" | PRIMARY_DISEASE == "Bipolar_Disorder_I" | PRIMARY_DISEASE == "Bipolar_Disorder_II" | PRIMARY_DISEASE == "Bipolar_Disorder_NOS")),
                      NA,
                      NA),
  "Controls" = c(nrow(df_merged %>% filter(PRIMARY_DISEASE == "Control")),
                 NA,
                 NA))

df_initial_summary_count <- rbindlist(list(df_initial_summary_count, list("Samples after initial filter", nrow(df_merged),
                                                                          nrow(df_merged %>% filter(PRIMARY_DISEASE == "Schizophrenia")),
                                                                          nrow(df_merged %>% filter(PRIMARY_DISEASE == "Bipolar_Disorder" | PRIMARY_DISEASE == "Bipolar_Disorder_I" | PRIMARY_DISEASE == "Bipolar_Disorder_II" | PRIMARY_DISEASE == "Bipolar_Disorder_NOS")),
                                                                          nrow(df_merged %>% filter(PRIMARY_DISEASE == "Control")))), use.names=FALSE)
fwrite(df_initial_summary_count, file='Analyses/initial_sample_count.tsv', quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')

df <- df_merged
names(df) <- gsub("sample_qc\\.", "", names(df))
names(df) <- gsub("phenotype\\.", "", names(df))

df_out <- filter(df, call_rate > T_sample_callRate) %>%
  filter(dp_stats.mean > T_dpMean) %>%
  filter(gq_stats.mean > T_gqMean)

df_out <- df_out %>% select(s)
print(dim(df_out))

fwrite(df_out, file=SAMPLE_LIST_INITIAL_QC, quote=FALSE, row.names=FALSE, col.names=FALSE)

# Create the table too
df_summary_count <- data.table(
  Filter = c("Samples after initial filter",
             paste0("Sample call rate < ", T_sample_callRate),
             paste0("Mean DP < ", T_dpMean),
             paste0("Mean GQ < ", T_gqMean),
             "Samples after sample QC filters"),
  Samples = c(nrow(df),
              nrow(filter(df, call_rate <= T_sample_callRate)),
              nrow(filter(df, dp_stats.mean <= T_dpMean)),
              nrow(filter(df, gq_stats.mean <= T_gqMean)),
              nrow(df_out)),
  "Schizophrenia Cases" = c(nrow(df %>% filter(PRIMARY_DISEASE == "Schizophrenia")),
                      nrow(filter(df, call_rate <= T_sample_callRate) %>% filter(PRIMARY_DISEASE == "Schizophrenia")),
                      nrow(filter(df, dp_stats.mean <= T_dpMean) %>% filter(PRIMARY_DISEASE == "Schizophrenia")),
                      nrow(filter(df, gq_stats.mean <= T_gqMean) %>% filter(PRIMARY_DISEASE == "Schizophrenia")),
                      length(which(df_out$s %in% (df %>% filter(PRIMARY_DISEASE == "Schizophrenia"))$s))),
  "Bipolar Cases" = c(nrow(df %>% filter(PRIMARY_DISEASE == "Bipolar_Disorder" | PRIMARY_DISEASE == "Bipolar_Disorder_I" | PRIMARY_DISEASE == "Bipolar_Disorder_II" | PRIMARY_DISEASE == "Bipolar_Disorder_NOS")),
                      nrow(filter(df, call_rate <= T_sample_callRate) %>% filter(PRIMARY_DISEASE == "Bipolar_Disorder" | PRIMARY_DISEASE == "Bipolar_Disorder_I" | PRIMARY_DISEASE == "Bipolar_Disorder_II" | PRIMARY_DISEASE == "Bipolar_Disorder_NOS")),
                      nrow(filter(df, dp_stats.mean <= T_dpMean) %>% filter(PRIMARY_DISEASE == "Bipolar_Disorder" | PRIMARY_DISEASE == "Bipolar_Disorder_I" | PRIMARY_DISEASE == "Bipolar_Disorder_II" | PRIMARY_DISEASE == "Bipolar_Disorder_NOS")),
                      nrow(filter(df, gq_stats.mean <= T_gqMean) %>% filter(PRIMARY_DISEASE == "Bipolar_Disorder" | PRIMARY_DISEASE == "Bipolar_Disorder_I" | PRIMARY_DISEASE == "Bipolar_Disorder_II" | PRIMARY_DISEASE == "Bipolar_Disorder_NOS")),
                      length(which(df_out$s %in% (df %>% filter(PRIMARY_DISEASE == "Bipolar_Disorder" | PRIMARY_DISEASE == "Bipolar_Disorder_I" | PRIMARY_DISEASE == "Bipolar_Disorder_II" | PRIMARY_DISEASE == "Bipolar_Disorder_NOS"))$s))),
  Controls = c(nrow(df %>% filter(PRIMARY_DISEASE == "Control")),
               nrow(filter(df, call_rate <= T_sample_callRate) %>% filter(PRIMARY_DISEASE == "Control")),
               nrow(filter(df, dp_stats.mean <= T_dpMean) %>% filter(PRIMARY_DISEASE == "Control")),
               nrow(filter(df, gq_stats.mean <= T_gqMean) %>% filter(PRIMARY_DISEASE == "Control")),
               length(which(df_out$s %in% (df %>% filter(PRIMARY_DISEASE == "Control"))$s))))

fwrite(df_summary_count, file='Analyses/sample_count.tsv', quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')
  