library(data.table)

INITIAL_VARIANT_QC_FILE <- 'Analyses/prefilter_metrics.tsv'
INITIAL_VARIANT_LIST <- 'Analyses/prefilter.keep.variant_list'

dt <- fread(INITIAL_VARIANT_QC_FILE)
n_initial_variants <- nrow(fread(INITIAL_VARIANT_LIST))

summary_df <- data.frame(
  Filter=c(
    "Variants with < 7 alleles",
    "Failing VQSR",
    "Invariant sites after initial variant and genotype filters",
    "Variants after initial filtering"
  ),
  Variants=c(
    nrow(dt),
    sum(dt$fail_VQSR),
    nrow(dt) - sum(dt$fail_VQSR) - n_initial_variants,
    n_initial_variants
  )
)

fwrite(summary_df, "Analyses/summary_variant_table.tsv", row.names=FALSE, sep='\t')
