library(data.table)

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

#association between URV variant types and schizophrenia 
glm_urv_ptv <- glm(outcome ~ n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv, family = binomial)
summary(glm_urv_ptv)
Sample <- "Total European"
Beta <- summary(glm_urv_ptv)$coefficients[2,1]
SE <- summary(glm_urv_ptv)$coefficients[2,2]
z_value <- summary(glm_urv_ptv)$coefficients[2,3]
p_value <- summary(glm_urv_ptv)$coefficients[2,4]
TotEUR_urv_ptv_est <- cbind.data.frame(Sample, Beta, SE, z_value, p_value)

glm_urv_syn <- glm(outcome ~ n_URV_synonymous + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv, family = binomial)
summary(glm_urv_syn)

#exome wide count as covariate
glm_urv_ptv <- glm(outcome ~ n_URV_PTV + n_URV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv, family = binomial)
summary(glm_urv_ptv)

glm_urv_syn <- glm(outcome ~ n_URV_synonymous + n_URV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv, family = binomial)
summary(glm_urv_syn)

#subset - include other variables and covariates 
mergeddata_urv_ptv <- subset(mergeddata_urv[,c("s","phenotype.SEX","outcome",colnames(mergeddata_urv)[grep("n_URV_PTV|pca.PC",colnames(mergeddata_urv))])])
mergeddata_urv_syn <- subset(mergeddata_urv[,c("s","phenotype.SEX","outcome",colnames(mergeddata_urv)[grep("n_URV_synonymous|pca.PC",colnames(mergeddata_urv))])])

#URV PTVs with exome-wide counts as covariates####
#run loop for all ptv variant burden scores across all genesets 
#outcome
out_start=3
out_end=3
out_nvar=out_end-out_start+1

out_variable=rep(NA, out_nvar)
out_beta=rep(NA, out_nvar)
out_se = rep(NA, out_nvar)
out_z = rep(NA, out_nvar)
out_OR = rep(NA, out_nvar)
out_p=rep(NA, out_nvar)

# exposure
exp_start=14
exp_end=125
exp_nvar=exp_end-exp_start+1

exp_variable=rep(NA, exp_nvar)
exp_beta=rep(NA, exp_nvar)
exp_se = rep(NA, out_nvar)
exp_z = rep(NA, out_nvar)
exp_OR = rep(NA, out_nvar)
exp_p=rep(NA, exp_nvar)

number=1

for (i in out_start:out_end) {
  outcome = colnames(mergeddata_urv_ptv)[i]
  for (j in exp_start:exp_end){
    exposure = colnames(mergeddata_urv_ptv)[j]
    x <- glm(outcome ~ get(exposure)+ n_URV_PTV + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_ptv, family = binomial)
    beta <- summary(x)$coefficients[,1]
    se <- summary(x)$coefficients[,2]
    z <- summary(x)$coefficients[,3]
    OR <- exp(beta)
    p <- summary(x)$coefficients[,4]
    
    out_beta[number] <- as.numeric(beta[2])
    out_se[number] <- as.numeric(se[2])
    out_z[number] <- as.numeric(z[2])
    out_OR[number] <- as.numeric(OR[2])
    out_p[number] <- as.numeric(p[2])
    out_variable[number] <- outcome
    number <- number + 1
    
    exp_beta[number] <- as.numeric(beta[2])
    exp_se[number] <- as.numeric(se[2])
    exp_z[number] <- as.numeric(z[2])
    exp_OR[number] <- as.numeric(OR[2])
    exp_p[number] <- as.numeric(p[2])
    exp_variable[number] <- exposure
    number <- number + 1
    
  }
}


outcome <- data.frame(out_variable, out_beta, out_se, out_z, out_OR, out_p)
exposure <- data.frame(exp_variable, exp_beta, exp_se, exp_z, exp_OR, exp_p)

#manage dataframes
library(tidyverse)

outcome <- outcome %>% 
  rename(
    variable = out_variable,
    beta = out_beta,
    se = out_se,
    z = out_z,
    OR = out_OR,
    p = out_p,
  )
exposure <- exposure %>% 
  rename(
    variable = exp_variable,
    beta = exp_beta,
    se = exp_se,
    z = exp_z,
    OR = exp_OR,
    p = exp_p,
  )
all <- rbind(outcome, exposure)
all <- na.omit(all)

#correct dataframe
dv <- all[c(1:112),]
iv <- all[c(113:224),]
dv$iv <- iv$variable

dv <- dv[c("variable", "iv", "beta", "se", "z", "OR", "p")]
dv$fdr_p <- p.adjust(dv$p, method = "fdr")
dv$significant <- ifelse(dv$fdr_p < 0.05, "yes", "no")
dv$outcome <- "Schizophrenia"
names(dv)[2] <- "geneset"
dv <- dv[c(10,2,3,4,5,6,7,8,9)]
table(dv$significant)

write.table(dv, file = "scz_urv_ptv_geneset_estimates_ewptvcov.tsv", sep="\t",dec = ".", row.names = F )

#URV synonymous variants  with exome-wide counts as covariates####
#run loop for all synonymous variant burden scores across all genesets 
#outcome
out_start=3
out_end=3
out_nvar=out_end-out_start+1

out_variable=rep(NA, out_nvar)
out_beta=rep(NA, out_nvar)
out_se = rep(NA, out_nvar)
out_z = rep(NA, out_nvar)
out_OR = rep(NA, out_nvar)
out_p=rep(NA, out_nvar)

# exposure
exp_start=14
exp_end=125
exp_nvar=exp_end-exp_start+1

exp_variable=rep(NA, exp_nvar)
exp_beta=rep(NA, exp_nvar)
exp_se = rep(NA, out_nvar)
exp_z = rep(NA, out_nvar)
exp_OR = rep(NA, out_nvar)
exp_p=rep(NA, exp_nvar)

number=1

for (i in out_start:out_end) {
  outcome = colnames(mergeddata_urv_syn)[i]
  for (j in exp_start:exp_end){
    exposure = colnames(mergeddata_urv_syn)[j]
    x <- glm(outcome ~ get(exposure)+ n_URV_synonymous + pca.PC1 + pca.PC2 + pca.PC3 + pca.PC4 + pca.PC5 + pca.PC6 + pca.PC7 + pca.PC8 + pca.PC9 + pca.PC10 + phenotype.SEX, data = mergeddata_urv_syn, family = binomial)
    beta <- summary(x)$coefficients[,1]
    se <- summary(x)$coefficients[,2]
    z <- summary(x)$coefficients[,3]
    OR <- exp(beta)
    p <- summary(x)$coefficients[,4]
    
    out_beta[number] <- as.numeric(beta[2])
    out_se[number] <- as.numeric(se[2])
    out_z[number] <- as.numeric(z[2])
    out_OR[number] <- as.numeric(OR[2])
    out_p[number] <- as.numeric(p[2])
    out_variable[number] <- outcome
    number <- number + 1
    
    exp_beta[number] <- as.numeric(beta[2])
    exp_se[number] <- as.numeric(se[2])
    exp_z[number] <- as.numeric(z[2])
    exp_OR[number] <- as.numeric(OR[2])
    exp_p[number] <- as.numeric(p[2])
    exp_variable[number] <- exposure
    number <- number + 1
    
  }
}


outcome <- data.frame(out_variable, out_beta, out_se, out_z, out_OR, out_p)
exposure <- data.frame(exp_variable, exp_beta, exp_se, exp_z, exp_OR, exp_p)

#manage dataframes
library(tidyverse)

outcome <- outcome %>% 
  rename(
    variable = out_variable,
    beta = out_beta,
    se = out_se,
    z = out_z,
    OR = out_OR,
    p = out_p,
  )
exposure <- exposure %>% 
  rename(
    variable = exp_variable,
    beta = exp_beta,
    se = exp_se,
    z = exp_z,
    OR = exp_OR,
    p = exp_p,
  )
all <- rbind(outcome, exposure)
all <- na.omit(all)

#correct dataframe
dv <- all[c(1:112),]
iv <- all[c(113:224),]
dv$iv <- iv$variable

dv <- dv[c("variable", "iv", "beta", "se", "z", "OR", "p")]
dv$fdr_p <- p.adjust(dv$p, method = "fdr")
dv$significant <- ifelse(dv$fdr_p < 0.05, "yes", "no")
dv$outcome <- "Schizophrenia"
names(dv)[2] <- "geneset"
dv <- dv[c(10,2,3,4,5,6,7,8,9)]
table(dv$significant)

write.table(dv, file = "scz_urv_syn_geneset_estimates_ewsyncov.tsv", sep="\t",dec = ".", row.names = F )