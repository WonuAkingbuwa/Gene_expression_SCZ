library(data.table)
library(openxlsx)

files <- list.files("/media/veracrypt10/Analyses/annotations/URVs_genesets/")
common <- fread(files[1], sep='\t', header=TRUE, data.table=FALSE, select = c(1,4:6,38,70:79))    # columns that are common to all files 
readfun <- function(fn) {
  x <- fread(fn, sep='\t', header=TRUE, data.table=FALSE, select = c(80:88))
  colnames(x) <- paste0(colnames(x), paste0("__",gsub('\\.', '_', fn)))
  return(x)
}
df <- do.call(cbind,lapply(files,readfun))
df <- cbind(common,df)

#check sex information is accurate 
table(df$imputesex.impute_sex.is_female)
table(df$phenotype.SEX)

#import overall URV file to use to for exclusions
#import overall URV counts
df_urv <- fread('/media/veracrypt10/Analyses/annotations/URVs_after_QC.tsv', sep='\t', header=TRUE, data.table=FALSE)
df_urv$n_URV <- df_urv$n_URV_SNP + df_urv$n_URV_indel
df_urv <- df_urv[c(1,80:89)]

#merge both files 
mergeddata_urv <- merge(df, df_urv, by='s')

#add dummy variable outcome column coded as 1 or 0
mergeddata_urv$outcome <- ifelse(mergeddata_urv$phenotype.ANALYSIS_CAT == "Case", 1, 0)

#filter to just people with schizophrenia 
mergeddata_urv <- subset(mergeddata_urv, mergeddata_urv$phenotype.PRIMARY_DISEASE == "Control" | mergeddata_urv$phenotype.PRIMARY_DISEASE == "Schizophrenia")

#exclude individuals with more than 4 median absolute deviations from the study and ethnicity-specific median number of ultra-rare variants
if (mad(mergeddata_urv$n_URV) == 0 )
{
  mergeddata_urv <- mergeddata_urv[(mergeddata_urv$n_URV < median(mergeddata_urv$n_URV) + 4*mad(mergeddata_urv$n_URV)) & (mergeddata_urv$n_URV > median(mergeddata_urv$n_URV) - 4*mad(mergeddata_urv$n_URV)),]
  
  print(paste0("N after exclusions:",nrow(mergeddata_urv)))
} else
{
  mergeddata_urv <- mergeddata_urv[(mergeddata_urv$n_URV < median(mergeddata_urv$n_URV) + 4*sd(mergeddata_urv$n_URV)) & (mergeddata_urv$n_URV > median(mergeddata_urv$n_URV) - 4*sd(mergeddata_urv$n_URV)),]
  
  print(paste0("N after exclusions:",nrow(mergeddata_urv)))
} 

#subset to just ptv/synonymous results - include other variables and covariates 
mergeddata_urv_ptv <- subset(mergeddata_urv[,c("s","phenotype.SEX","outcome",colnames(mergeddata_urv)[grep("n_URV_PTV|pca.PC",colnames(mergeddata_urv))])])
mergeddata_urv_syn <- subset(mergeddata_urv[,c("s","phenotype.SEX","outcome",colnames(mergeddata_urv)[grep("n_URV_synonymous|pca.PC",colnames(mergeddata_urv))])])

#to match the conditional analyses in MAGMA for the common variants, we will run regressions for the PIx gene sets that were significant in previous 
#analyses. first we will include URV counts from their corresponding brain cell gene sets as covariates, then PI genes, then both. 

####2) conditional analyses; exome-wide correction with n_URV_PTV####
#2.1) corresponding brain cells as covariates ####
#mouse brain cells
#PI x striatal interneurons
pistriint_striint_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxStriatal_Interneuron_intervals_bed_tsv + n_URV_PTV__URVs_Striatal_Interneuron_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(pistriint_striint_cov)
Category <- "Mouse brain cells"
Geneset <- "PI x Striatal Interneurons"
Beta <- summary(pistriint_striint_cov)$coefficients[2,1]
SE <- summary(pistriint_striint_cov)$coefficients[2,2]
z_value <- summary(pistriint_striint_cov)$coefficients[2,3]
p_value <- summary(pistriint_striint_cov)$coefficients[2,4]
results_pistriint_striint_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#PI x pyramidal S1
pipyss_pyss_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxpyramidal_SS_intervals_bed_tsv + n_URV_PTV__URVs_pyramidal_SS_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(pipyss_pyss_cov)
Category <- "Mouse brain cells"
Geneset <- "PI x Pyramidal (S1)"
Beta <- summary(pipyss_pyss_cov)$coefficients[2,1]
SE <- summary(pipyss_pyss_cov)$coefficients[2,2]
z_value <- summary(pipyss_pyss_cov)$coefficients[2,3]
p_value <- summary(pipyss_pyss_cov)$coefficients[2,4]
results_pipyss_pyss_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#PI x pyramidal CA1
pipyca1_pyca1_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxpyramidal_CA1_intervals_bed_tsv + n_URV_PTV__URVs_pyramidal_CA1_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(pipyca1_pyca1_cov)
Category <- "Mouse brain cells"
Geneset <- "PI x Pyramidal (CA1)"
Beta <- summary(pipyca1_pyca1_cov)$coefficients[2,1]
SE <- summary(pipyca1_pyca1_cov)$coefficients[2,2]
z_value <- summary(pipyca1_pyca1_cov)$coefficients[2,3]
p_value <- summary(pipyca1_pyca1_cov)$coefficients[2,4]
results_pipyca1_pyca1_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#PI x medium spiny neurons
pimsn_msn_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxMedium_Spiny_Neuron_intervals_bed_tsv + n_URV_PTV__URVs_Medium_Spiny_Neuron_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(pimsn_msn_cov)
Category <- "Mouse brain cells"
Geneset <- "PI x Medium Spiny Neuron"
Beta <- summary(pimsn_msn_cov)$coefficients[2,1]
SE <- summary(pimsn_msn_cov)$coefficients[2,2]
z_value <- summary(pimsn_msn_cov)$coefficients[2,3]
p_value <- summary(pimsn_msn_cov)$coefficients[2,4]
results_pimsn_msn_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#PI x interneurons
piint_int_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxinterneurons_intervals_bed_tsv + n_URV_PTV__URVs_interneurons_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(piint_int_cov)
Category <- "Mouse brain cells"
Geneset <- "PI x Interneurons"
Beta <- summary(piint_int_cov)$coefficients[2,1]
SE <- summary(piint_int_cov)$coefficients[2,2]
z_value <- summary(piint_int_cov)$coefficients[2,3]
p_value <- summary(piint_int_cov)$coefficients[2,4]
results_piint_int_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#PI x hypothalamic dopaminergic neurons
pihdn_hdn_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxHypothalamic_Dopaminergic_intervals_bed_tsv + n_URV_PTV__URVs_Hypothalamic_Dopaminergic_Neurons_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(pihdn_hdn_cov)
Category <- "Mouse brain cells"
Geneset <- "PI x Hypothalamic Dopaminergic Neurons"
Beta <- summary(pihdn_hdn_cov)$coefficients[2,1]
SE <- summary(pihdn_hdn_cov)$coefficients[2,2]
z_value <- summary(pihdn_hdn_cov)$coefficients[2,3]
p_value <- summary(pihdn_hdn_cov)$coefficients[2,4]
results_pihdn_hdn_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#PI x dopaminergic neuroblast
pidn_dn_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxDopaminergic_Neuroblast_intervals_bed_tsv + n_URV_PTV__URVs_Dopaminergic_Neuroblast_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(pidn_dn_cov)
Category <- "Mouse brain cells"
Geneset <- "PI x Dopaminergic Neuroblast"
Beta <- summary(pidn_dn_cov)$coefficients[2,1]
SE <- summary(pidn_dn_cov)$coefficients[2,2]
z_value <- summary(pidn_dn_cov)$coefficients[2,3]
p_value <- summary(pidn_dn_cov)$coefficients[2,4]
results_pidn_dn_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#human brain cells
#GABA2
pigaba2_gaba2_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxGABA2_intervals_bed_tsv + n_URV_PTV__URVs_GABA2_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(pigaba2_gaba2_cov)
Category <- "Human brain cells"
Geneset <- "PI x GABA2"
Beta <- summary(pigaba2_gaba2_cov)$coefficients[2,1]
SE <- summary(pigaba2_gaba2_cov)$coefficients[2,2]
z_value <- summary(pigaba2_gaba2_cov)$coefficients[2,3]
p_value <- summary(pigaba2_gaba2_cov)$coefficients[2,4]
results_pigaba2_gaba2_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#GABA1
pigaba1_gaba1_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxGABA1_intervals_bed_tsv + n_URV_PTV__URVs_GABA1_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(pigaba1_gaba1_cov)
Category <- "Human brain cells"
Geneset <- "PI x GABA1"
Beta <- summary(pigaba1_gaba1_cov)$coefficients[2,1]
SE <- summary(pigaba1_gaba1_cov)$coefficients[2,2]
z_value <- summary(pigaba1_gaba1_cov)$coefficients[2,3]
p_value <- summary(pigaba1_gaba1_cov)$coefficients[2,4]
results_pigaba1_gaba1_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#exPFC1
piexpfc1_expfc1_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxexPFC1_intervals_bed_tsv + n_URV_PTV__URVs_exPFC1_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(piexpfc1_expfc1_cov)
Category <- "Human brain cells"
Geneset <- "PI x exPFC1"
Beta <- summary(piexpfc1_expfc1_cov)$coefficients[2,1]
SE <- summary(piexpfc1_expfc1_cov)$coefficients[2,2]
z_value <- summary(piexpfc1_expfc1_cov)$coefficients[2,3]
p_value <- summary(piexpfc1_expfc1_cov)$coefficients[2,4]
results_piexpfc1_expfc1_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#exDG
piexdg_exdg_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxexDG_intervals_bed_tsv + n_URV_PTV__URVs_exDG_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(piexdg_exdg_cov)
Category <- "Human brain cells"
Geneset <- "PI x exDG"
Beta <- summary(piexdg_exdg_cov)$coefficients[2,1]
SE <- summary(piexdg_exdg_cov)$coefficients[2,2]
z_value <- summary(piexdg_exdg_cov)$coefficients[2,3]
p_value <- summary(piexdg_exdg_cov)$coefficients[2,4]
results_piexdg_exdg_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#exCA1
piexca1_exca1_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxexCA1_intervals_bed_tsv + n_URV_PTV__URVs_exCA1_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(piexca1_exca1_cov)
Category <- "Human brain cells"
Geneset <- "PI x exCA1"
Beta <- summary(piexca1_exca1_cov)$coefficients[2,1]
SE <- summary(piexca1_exca1_cov)$coefficients[2,2]
z_value <- summary(piexca1_exca1_cov)$coefficients[2,3]
p_value <- summary(piexca1_exca1_cov)$coefficients[2,4]
results_piexca1_exca1_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#ASC1
piasc1_asc1_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxASC1_intervals_bed_tsv + n_URV_PTV__URVs_ASC1_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(piasc1_asc1_cov)
Category <- "Human brain cells"
Geneset <- "PI x ASC1"
Beta <- summary(piasc1_asc1_cov)$coefficients[2,1]
SE <- summary(piasc1_asc1_cov)$coefficients[2,2]
z_value <- summary(piasc1_asc1_cov)$coefficients[2,3]
p_value <- summary(piasc1_asc1_cov)$coefficients[2,4]
results_piasc1_asc1_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#row bind and write out results table
results_corresponding_sets_cov <- rbind(results_pistriint_striint_cov, results_pipyss_pyss_cov, results_pipyca1_pyca1_cov, results_pimsn_msn_cov, results_piint_int_cov, results_pihdn_hdn_cov, results_pidn_dn_cov, 
                                        results_pigaba2_gaba2_cov, results_pigaba1_gaba1_cov, results_piexpfc1_expfc1_cov, results_piexdg_exdg_cov, results_piexca1_exca1_cov, results_piasc1_asc1_cov)
write.xlsx(results_corresponding_sets_cov, "/media/veracrypt10/Analyses/annotations/regression_results/conditional_results_corresponding_sets_cov.xlsx", row.names = F, col.name = T)

#2.2) PI genes as covariates ####
#mouse brain cells
#PI x striatal interneurons
pistriint_pi_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxStriatal_Interneuron_intervals_bed_tsv + n_URV_PTV__URVs_PI_genes_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(pistriint_pi_cov)
Category <- "Mouse brain cells"
Geneset <- "PI x Striatal Interneurons"
Beta <- summary(pistriint_pi_cov)$coefficients[2,1]
SE <- summary(pistriint_pi_cov)$coefficients[2,2]
z_value <- summary(pistriint_pi_cov)$coefficients[2,3]
p_value <- summary(pistriint_pi_cov)$coefficients[2,4]
results_pistriint_pi_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#PI x pyramidal S1
pipyss_pi_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxpyramidal_SS_intervals_bed_tsv + n_URV_PTV__URVs_PI_genes_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(pipyss_pi_cov)
Category <- "Mouse brain cells"
Geneset <- "PI x Pyramidal (S1)"
Beta <- summary(pipyss_pi_cov)$coefficients[2,1]
SE <- summary(pipyss_pi_cov)$coefficients[2,2]
z_value <- summary(pipyss_pi_cov)$coefficients[2,3]
p_value <- summary(pipyss_pi_cov)$coefficients[2,4]
results_pipyss_pi_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#PI x pyramidal CA1
pipyca1_pi_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxpyramidal_CA1_intervals_bed_tsv + n_URV_PTV__URVs_PI_genes_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(pipyca1_pi_cov)
Category <- "Mouse brain cells"
Geneset <- "PI x Pyramidal (CA1)"
Beta <- summary(pipyca1_pi_cov)$coefficients[2,1]
SE <- summary(pipyca1_pi_cov)$coefficients[2,2]
z_value <- summary(pipyca1_pi_cov)$coefficients[2,3]
p_value <- summary(pipyca1_pi_cov)$coefficients[2,4]
results_pipyca1_pi_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#PI x medium spiny neurons
pimsn_pi_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxMedium_Spiny_Neuron_intervals_bed_tsv + n_URV_PTV__URVs_PI_genes_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(pimsn_pi_cov)
Category <- "Mouse brain cells"
Geneset <- "PI x Medium Spiny Neuron"
Beta <- summary(pimsn_pi_cov)$coefficients[2,1]
SE <- summary(pimsn_pi_cov)$coefficients[2,2]
z_value <- summary(pimsn_pi_cov)$coefficients[2,3]
p_value <- summary(pimsn_pi_cov)$coefficients[2,4]
results_pimsn_pi_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#PI x interneurons
piint_pi_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxinterneurons_intervals_bed_tsv + n_URV_PTV__URVs_PI_genes_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(piint_pi_cov)
Category <- "Mouse brain cells"
Geneset <- "PI x Interneurons"
Beta <- summary(piint_pi_cov)$coefficients[2,1]
SE <- summary(piint_pi_cov)$coefficients[2,2]
z_value <- summary(piint_pi_cov)$coefficients[2,3]
p_value <- summary(piint_pi_cov)$coefficients[2,4]
results_piint_pi_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#PI x hypothalamic dopaminergic neurons
pihdn_pi_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxHypothalamic_Dopaminergic_intervals_bed_tsv + n_URV_PTV__URVs_PI_genes_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(pihdn_pi_cov)
Category <- "Mouse brain cells"
Geneset <- "PI x Hypothalamic Dopaminergic Neurons"
Beta <- summary(pihdn_pi_cov)$coefficients[2,1]
SE <- summary(pihdn_pi_cov)$coefficients[2,2]
z_value <- summary(pihdn_pi_cov)$coefficients[2,3]
p_value <- summary(pihdn_pi_cov)$coefficients[2,4]
results_pihdn_pi_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#PI x dopaminergic neuroblast
pidn_pi_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxDopaminergic_Neuroblast_intervals_bed_tsv + n_URV_PTV__URVs_PI_genes_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(pidn_pi_cov)
Category <- "Mouse brain cells"
Geneset <- "PI x Dopaminergic Neuroblast"
Beta <- summary(pidn_pi_cov)$coefficients[2,1]
SE <- summary(pidn_pi_cov)$coefficients[2,2]
z_value <- summary(pidn_pi_cov)$coefficients[2,3]
p_value <- summary(pidn_pi_cov)$coefficients[2,4]
results_pidn_pi_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#human brain cells
#GABA2
pigaba2_pi_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxGABA2_intervals_bed_tsv + n_URV_PTV__URVs_PI_genes_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(pigaba2_pi_cov)
Category <- "Human brain cells"
Geneset <- "PI x GABA2"
Beta <- summary(pigaba2_pi_cov)$coefficients[2,1]
SE <- summary(pigaba2_pi_cov)$coefficients[2,2]
z_value <- summary(pigaba2_pi_cov)$coefficients[2,3]
p_value <- summary(pigaba2_pi_cov)$coefficients[2,4]
results_pigaba2_pi_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#GABA1
pigaba1_pi_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxGABA1_intervals_bed_tsv + n_URV_PTV__URVs_PI_genes_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(pigaba1_pi_cov)
Category <- "Human brain cells"
Geneset <- "PI x GABA1"
Beta <- summary(pigaba1_pi_cov)$coefficients[2,1]
SE <- summary(pigaba1_pi_cov)$coefficients[2,2]
z_value <- summary(pigaba1_pi_cov)$coefficients[2,3]
p_value <- summary(pigaba1_pi_cov)$coefficients[2,4]
results_pigaba1_pi_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#exPFC1
piexpfc1_pi_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxexPFC1_intervals_bed_tsv + n_URV_PTV__URVs_PI_genes_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(piexpfc1_pi_cov)
Category <- "Human brain cells"
Geneset <- "PI x exPFC1"
Beta <- summary(piexpfc1_pi_cov)$coefficients[2,1]
SE <- summary(piexpfc1_pi_cov)$coefficients[2,2]
z_value <- summary(piexpfc1_pi_cov)$coefficients[2,3]
p_value <- summary(piexpfc1_pi_cov)$coefficients[2,4]
results_piexpfc1_pi_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#exDG
piexdg_pi_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxexDG_intervals_bed_tsv + n_URV_PTV__URVs_PI_genes_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(piexdg_pi_cov)
Category <- "Human brain cells"
Geneset <- "PI x exDG"
Beta <- summary(piexdg_pi_cov)$coefficients[2,1]
SE <- summary(piexdg_pi_cov)$coefficients[2,2]
z_value <- summary(piexdg_pi_cov)$coefficients[2,3]
p_value <- summary(piexdg_pi_cov)$coefficients[2,4]
results_piexdg_pi_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#exCA1
piexca1_pi_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxexCA1_intervals_bed_tsv + n_URV_PTV__URVs_PI_genes_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(piexca1_pi_cov)
Category <- "Human brain cells"
Geneset <- "PI x exCA1"
Beta <- summary(piexca1_pi_cov)$coefficients[2,1]
SE <- summary(piexca1_pi_cov)$coefficients[2,2]
z_value <- summary(piexca1_pi_cov)$coefficients[2,3]
p_value <- summary(piexca1_pi_cov)$coefficients[2,4]
results_piexca1_pi_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#ASC1
piasc1_pi_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxASC1_intervals_bed_tsv + n_URV_PTV__URVs_PI_genes_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(piasc1_pi_cov)
Category <- "Human brain cells"
Geneset <- "PI x ASC1"
Beta <- summary(piasc1_pi_cov)$coefficients[2,1]
SE <- summary(piasc1_pi_cov)$coefficients[2,2]
z_value <- summary(piasc1_pi_cov)$coefficients[2,3]
p_value <- summary(piasc1_pi_cov)$coefficients[2,4]
results_piasc1_pi_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#row bind and write out results table
results_pi_cov <- rbind(results_pistriint_pi_cov, results_pipyss_pi_cov, results_pipyca1_pi_cov, results_pimsn_pi_cov, results_piint_pi_cov, results_pihdn_pi_cov, results_pidn_pi_cov, 
                        results_pigaba2_pi_cov, results_pigaba1_pi_cov, results_piexpfc1_pi_cov, results_piexdg_pi_cov, results_piexca1_pi_cov, results_piasc1_pi_cov)
write.xlsx(results_pi_cov, "/media/veracrypt10/Analyses/annotations/regression_results/conditional_results_pi_cov.xlsx", row.names = F, col.name = T)

#2.3) PI genes + corresponding sets as covariates ####
#mouse brain cells
#PI x striatal interneurons
pistriint_both_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxStriatal_Interneuron_intervals_bed_tsv + n_URV_PTV__URVs_Striatal_Interneuron_intervals_bed_tsv + n_URV_PTV__URVs_PI_genes_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(pistriint_both_cov)
Category <- "Mouse brain cells"
Geneset <- "PI x Striatal Interneurons"
Beta <- summary(pistriint_both_cov)$coefficients[2,1]
SE <- summary(pistriint_both_cov)$coefficients[2,2]
z_value <- summary(pistriint_both_cov)$coefficients[2,3]
p_value <- summary(pistriint_both_cov)$coefficients[2,4]
results_pistriint_both_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#PI x pyramidal S1
pipyss_both_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxpyramidal_SS_intervals_bed_tsv + n_URV_PTV__URVs_pyramidal_SS_intervals_bed_tsv + n_URV_PTV__URVs_PI_genes_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(pipyss_both_cov)
Category <- "Mouse brain cells"
Geneset <- "PI x Pyramidal (S1)"
Beta <- summary(pipyss_both_cov)$coefficients[2,1]
SE <- summary(pipyss_both_cov)$coefficients[2,2]
z_value <- summary(pipyss_both_cov)$coefficients[2,3]
p_value <- summary(pipyss_both_cov)$coefficients[2,4]
results_pipyss_both_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#PI x pyramidal CA1
pipyca1_both_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxpyramidal_CA1_intervals_bed_tsv + n_URV_PTV__URVs_pyramidal_CA1_intervals_bed_tsv + n_URV_PTV__URVs_PI_genes_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(pipyca1_both_cov)
Category <- "Mouse brain cells"
Geneset <- "PI x Pyramidal (CA1)"
Beta <- summary(pipyca1_both_cov)$coefficients[2,1]
SE <- summary(pipyca1_both_cov)$coefficients[2,2]
z_value <- summary(pipyca1_both_cov)$coefficients[2,3]
p_value <- summary(pipyca1_both_cov)$coefficients[2,4]
results_pipyca1_both_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#PI x medium spiny neurons
pimsn_both_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxMedium_Spiny_Neuron_intervals_bed_tsv + n_URV_PTV__URVs_Medium_Spiny_Neuron_intervals_bed_tsv + n_URV_PTV__URVs_PI_genes_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(pimsn_both_cov)
Category <- "Mouse brain cells"
Geneset <- "PI x Medium Spiny Neuron"
Beta <- summary(pimsn_both_cov)$coefficients[2,1]
SE <- summary(pimsn_both_cov)$coefficients[2,2]
z_value <- summary(pimsn_both_cov)$coefficients[2,3]
p_value <- summary(pimsn_both_cov)$coefficients[2,4]
results_pimsn_both_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#PI x interneurons
piint_both_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxinterneurons_intervals_bed_tsv + n_URV_PTV__URVs_interneurons_intervals_bed_tsv + n_URV_PTV__URVs_PI_genes_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(piint_both_cov)
Category <- "Mouse brain cells"
Geneset <- "PI x Interneurons"
Beta <- summary(piint_both_cov)$coefficients[2,1]
SE <- summary(piint_both_cov)$coefficients[2,2]
z_value <- summary(piint_both_cov)$coefficients[2,3]
p_value <- summary(piint_both_cov)$coefficients[2,4]
results_piint_both_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#PI x hypothalamic dopaminergic neurons
pihdn_both_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxHypothalamic_Dopaminergic_intervals_bed_tsv + n_URV_PTV__URVs_Hypothalamic_Dopaminergic_Neurons_intervals_bed_tsv + n_URV_PTV__URVs_PI_genes_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(pihdn_both_cov)
Category <- "Mouse brain cells"
Geneset <- "PI x Hypothalamic Dopaminergic Neurons"
Beta <- summary(pihdn_both_cov)$coefficients[2,1]
SE <- summary(pihdn_both_cov)$coefficients[2,2]
z_value <- summary(pihdn_both_cov)$coefficients[2,3]
p_value <- summary(pihdn_both_cov)$coefficients[2,4]
results_pihdn_both_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#PI x dopaminergic neuroblast
pidn_both_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxDopaminergic_Neuroblast_intervals_bed_tsv + n_URV_PTV__URVs_Dopaminergic_Neuroblast_intervals_bed_tsv + n_URV_PTV__URVs_PI_genes_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(pidn_both_cov)
Category <- "Mouse brain cells"
Geneset <- "PI x Dopaminergic Neuroblast"
Beta <- summary(pidn_both_cov)$coefficients[2,1]
SE <- summary(pidn_both_cov)$coefficients[2,2]
z_value <- summary(pidn_both_cov)$coefficients[2,3]
p_value <- summary(pidn_both_cov)$coefficients[2,4]
results_pidn_both_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#human brain cells
#GABA2
pigaba2_both_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxGABA2_intervals_bed_tsv + n_URV_PTV__URVs_GABA2_intervals_bed_tsv + n_URV_PTV__URVs_PI_genes_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(pigaba2_both_cov)
Category <- "Human brain cells"
Geneset <- "PI x GABA2"
Beta <- summary(pigaba2_both_cov)$coefficients[2,1]
SE <- summary(pigaba2_both_cov)$coefficients[2,2]
z_value <- summary(pigaba2_both_cov)$coefficients[2,3]
p_value <- summary(pigaba2_both_cov)$coefficients[2,4]
results_pigaba2_both_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#GABA1
pigaba1_both_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxGABA1_intervals_bed_tsv + n_URV_PTV__URVs_GABA1_intervals_bed_tsv + n_URV_PTV__URVs_PI_genes_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(pigaba1_both_cov)
Category <- "Human brain cells"
Geneset <- "PI x GABA1"
Beta <- summary(pigaba1_both_cov)$coefficients[2,1]
SE <- summary(pigaba1_both_cov)$coefficients[2,2]
z_value <- summary(pigaba1_both_cov)$coefficients[2,3]
p_value <- summary(pigaba1_both_cov)$coefficients[2,4]
results_pigaba1_both_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#exPFC1
piexpfc1_both_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxexPFC1_intervals_bed_tsv + n_URV_PTV__URVs_exPFC1_intervals_bed_tsv + n_URV_PTV__URVs_PI_genes_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(piexpfc1_both_cov)
Category <- "Human brain cells"
Geneset <- "PI x exPFC1"
Beta <- summary(piexpfc1_both_cov)$coefficients[2,1]
SE <- summary(piexpfc1_both_cov)$coefficients[2,2]
z_value <- summary(piexpfc1_both_cov)$coefficients[2,3]
p_value <- summary(piexpfc1_both_cov)$coefficients[2,4]
results_piexpfc1_both_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#exDG
piexdg_both_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxexDG_intervals_bed_tsv + n_URV_PTV__URVs_exDG_intervals_bed_tsv + n_URV_PTV__URVs_PI_genes_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(piexdg_both_cov)
Category <- "Human brain cells"
Geneset <- "PI x exDG"
Beta <- summary(piexdg_both_cov)$coefficients[2,1]
SE <- summary(piexdg_both_cov)$coefficients[2,2]
z_value <- summary(piexdg_both_cov)$coefficients[2,3]
p_value <- summary(piexdg_both_cov)$coefficients[2,4]
results_piexdg_both_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#exCA1
piexca1_both_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxexCA1_intervals_bed_tsv + n_URV_PTV__URVs_exCA1_intervals_bed_tsv + n_URV_PTV__URVs_PI_genes_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(piexca1_both_cov)
Category <- "Human brain cells"
Geneset <- "PI x exCA1"
Beta <- summary(piexca1_both_cov)$coefficients[2,1]
SE <- summary(piexca1_both_cov)$coefficients[2,2]
z_value <- summary(piexca1_both_cov)$coefficients[2,3]
p_value <- summary(piexca1_both_cov)$coefficients[2,4]
results_piexca1_both_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#ASC1
piasc1_both_cov <- glm(outcome ~ n_URV_PTV__URVs_PIxASC1_intervals_bed_tsv + n_URV_PTV__URVs_ASC1_intervals_bed_tsv + n_URV_PTV__URVs_PI_genes_intervals_bed_tsv + n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
summary(piasc1_both_cov)
Category <- "Human brain cells"
Geneset <- "PI x ASC1"
Beta <- summary(piasc1_both_cov)$coefficients[2,1]
SE <- summary(piasc1_both_cov)$coefficients[2,2]
z_value <- summary(piasc1_both_cov)$coefficients[2,3]
p_value <- summary(piasc1_both_cov)$coefficients[2,4]
results_piasc1_both_cov <- cbind.data.frame(Category, Geneset, Beta, SE, z_value, p_value)

#row bind and write out results table
results_both_cov <- rbind(results_pistriint_both_cov, results_pipyss_both_cov, results_pipyca1_both_cov, results_pimsn_both_cov, results_piint_both_cov, results_pihdn_both_cov, results_pidn_both_cov, 
                          results_pigaba2_both_cov, results_pigaba1_both_cov, results_piexpfc1_both_cov, results_piexdg_both_cov, results_piexca1_both_cov, results_piasc1_both_cov)
write.xlsx(results_both_cov, "/media/veracrypt10/Analyses/annotations/regression_results/conditional_results_both_cov.xlsx", row.names = F, col.name = T)
























