library(ggplot2)
library(dplyr)
library(ggsci)
library(data.table)
library(randomForest)

source("Scripts/r_options_Wonu.R")
source('Scripts/pretty_plotting.R')
source("Scripts/helpful_functions.R")

# Throughout, we only consider the strictly defined European set.
save_figures <- TRUE

URV_FILE <- 'annotations/URVs_after_QC.tsv'

df_URVs <- fread(URV_FILE, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE)
names(df_URVs)[names(df_URVs) == 'Samples'] <- 's'

# All the European samples - same as strictly defined Europeans
df_1kg_EUR <- fread(PCA_EUR_1KG_SCORES, sep='\t', header=TRUE, data.table=FALSE)
df_1kg_EUR <- df_1kg_EUR[sample(nrow(df_1kg_EUR), replace=FALSE),]
df_1kg_EUR <- merge(df_1kg_EUR, df_URVs, all.x=TRUE)

#plot PC1 and PC2 to exclude Finnish outliers#### 
plot(df_1kg_EUR$PC1, df_1kg_EUR$PC2, col=as.factor(df_1kg_EUR$POPULATION))
legend('topright', legend = levels(as.factor(df_1kg_EUR$POPULATION)), col = 1:7, cex = 0.6, pch = 1, pt.cex = 2)

ellipse.1 <- ellipse(2.5, 0.125, -0.03, 0.085, 0.09)
lines(ellipse.1[,1], ellipse.1[,2], lwd=2, col='indianred3')

in_ell_1 <- in_ellipse(df_1kg_EUR$PC1, df_1kg_EUR$PC2, 2.5, 0.125, -0.03, 0.085, 0.09)
sum(in_ell_1 == TRUE) #will include 1KG individuals as well

# Create labellings
FIN <- 2 * in_ell_1
table(FIN) 

df_1kg_EUR$FIN <- factor(FIN, labels=c("not FIN", "FIN"))
table(df_1kg_EUR$FIN) 

#Next, what we think are Northern Swedish outliers#### 
ellipse.2 <- ellipse(2.5, 0.03, 0.09, 0.05, 0.07)
lines(ellipse.2[,1], ellipse.2[,2], lwd=2, col='indianred3')

in_ell_2 <- in_ellipse(df_1kg_EUR$PC1, df_1kg_EUR$PC2, 2.5, 0.03, 0.09, 0.05, 0.07)
sum(in_ell_2 == TRUE) 

# Create labellings
NorSwe <- 2 * in_ell_2
table(NorSwe) 

df_1kg_EUR$NorSwe <- factor(NorSwe, labels=c("not NorSwe", "NorSwe"))
table(df_1kg_EUR$NorSwe) 

#other potential EUR outliers####
ellipse.3 <- ellipse(2.5, -0.07, -0.06, 0.04, 0.03)
lines(ellipse.3[,1], ellipse.3[,2], lwd=2, col='indianred3')

in_ell_3 <- in_ellipse(df_1kg_EUR$PC1, df_1kg_EUR$PC2, 2.5, -0.07, -0.06, 0.04, 0.03)
sum(in_ell_3 == TRUE) 

# Create labellings
OthEx <- 2 * in_ell_3
table(OthEx) #0=12313, 2=264

df_1kg_EUR$OthEx <- factor(OthEx, labels=c("not OthEx", "OthEx"))
table(df_1kg_EUR$OthEx) # not OthEx=12313, OthEx=264

#check exclusions in just my data####
#subset to just my data
df_EUR <- subset(df_1kg_EUR, df_1kg_EUR$POPULATION == "Case" | df_1kg_EUR$POPULATION == "Control")

table(df_EUR$FIN) 
table(df_EUR$NorSwe) 
table(df_EUR$OthEx) 

#make column with total outliers
df_EUR$outlier <- ifelse(df_EUR$FIN == "FIN" | df_EUR$NorSwe == "NorSwe" | df_EUR$OthEx == "OthEx", "outlier", "not outlier")
table(df_EUR$outlier) 

#make exclusionary list of each of the three groups + single list of all three
finnish <- df_EUR %>% filter(df_EUR$FIN == "FIN")
finnish <- finnish %>% select(s)
fwrite(finnish, file="Analyses/finnish_pca.remove_sample_list", quote=FALSE, row.names=FALSE, col.names=FALSE)

swedish <- df_EUR %>% filter(df_EUR$NorSwe == "NorSwe")
swedish <- swedish %>% select(s)
fwrite(swedish, file="Analyses/NorSwe_pca.remove_sample_list", quote=FALSE, row.names=FALSE, col.names=FALSE)

other <- df_EUR %>% filter(df_EUR$OthEx == "OthEx")
other <- other %>% select(s)
fwrite(other, file="Analyses/otherEUR_pca.remove_sample_list", quote=FALSE, row.names=FALSE, col.names=FALSE)

all_ouliers <- df_EUR %>% filter(df_EUR$outlier == "outlier")
all_ouliers <- all_ouliers %>% select(s)
fwrite(all_ouliers, file="Analyses/EURoutliers_pca.remove_sample_list", quote=FALSE, row.names=FALSE, col.names=FALSE)




