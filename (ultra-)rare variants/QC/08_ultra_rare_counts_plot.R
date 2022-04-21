
library(hexbin)
library(gridExtra)
library(ggplot2)
library(ggExtra)
library(data.table)

# Load plotting functions:
source('Scripts/pretty_plotting.R')

# Get file locations and plotting locations.
source("Scripts/r_options_Wonu.R")

URV_FILE <- 'annotations/URVs.tsv'

df <- fread(URV_FILE, sep='\t', header=TRUE, data.table=FALSE)
df <- df[sample(nrow(df), replace=FALSE),]

# Scatters of URV-SNPs against URV-Indels.
d <- ggplot(df, aes(x=n_URV_SNP, y=n_URV_indel, colour=phenotype.SAMP_SOURCE)) + 
  geom_point(size=0.5, alpha=0.5) + 
  scale_color_d3('category10') + theme_minimal() + labs(x='n URV SNPs', y='n URV indels', color='Batch')
d <- ggExtra::ggMarginal(d, type = "density",
                    xparams = list(adjust=1), yparams=list(adjust=1.5))
print(d)

save_figures <- TRUE

y_labels <- c('Batch', 'Location')
y_label_batch <- c('', '')
titles <- c('Number of Singletons split by Batch and coloured by Location',
	'Number of Singletons split by Location and coloured by Phenotype')
titles <- c('', '')

create_pretty_boxplots(df, aes(x=phenotype.SAMP_SOURCE, y=n_URV_SNP), aes(colour=factor(phenotype.SITE)), 
	y_label=y_label_batch[1], x_label='Number of Singletons', key_label='Location',
	title=titles[1], legend=TRUE, save_figure=save_figures,  file=paste0(PLOTS,'08_URVs_by_batch'), xlim=c(0,700),
	threshold=T_nURVSNP)
create_pretty_boxplots(df, aes(x=phenotype.SITE, y=n_URV_SNP), aes(colour=factor(phenotype.PRIMARY_DISEASE)),
	y_label=y_label_batch[2], x_label='Number of Singletons', key_label='Phenotype',
	title=titles[2], legend=TRUE, save_figure=save_figures,  file=paste0(PLOTS,'08_URVs_by_location'),  xlim=c(0,700),
	threshold=T_nURVSNP)

# Should remove individuals with more than 300 singletons.

create_pretty_boxplots(df, aes(x=phenotype.SAMP_SOURCE, y=n_URV_indel), aes(colour=factor(phenotype.SITE)), 
	y_label=y_label_batch[1], x_label='Number of Singletons', key_label='Location',
	title=titles[1], legend=TRUE, save_figure=save_figures,  file=paste0(PLOTS,'08_URV_indels_by_batch'), xlim=c(0,100),
	threshold=T_nURVIndel)
create_pretty_boxplots(df, aes(x=phenotype.SITE, y=n_URV_indel), aes(colour=factor(phenotype.PRIMARY_DISEASE)),
	y_label=y_label_batch[2], x_label='Number of Singletons', key_label='Phenotype',
	title=titles[2], legend=TRUE, save_figure=save_figures,  file=paste0(PLOTS,'08_URV_indels_by_location'),  xlim=c(0,100),
	threshold=T_nURVIndel)
