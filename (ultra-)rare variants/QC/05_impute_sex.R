rm(list=ls())
library(ggplot2)
library(ggsci)
library(dplyr)
library(data.table)

# File locations and plotting locations defined in r_options_Wonu.R
source("r_options_Wonu.R")
source("pretty_plotting.R")


IMPUTESEX_FILE <- "Analyses/imputesex.tsv"
Y_NCALLED_FILE <- "Analyses/ycalled.tsv"

SEXCHECK_LIST <- 'Analyses/sexcheck.remove.sample_list'

df <- fread(IMPUTESEX_FILE, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE) %>%
  mutate(imputed_sex=as.factor(ifelse(impute_sex.is_female == TRUE, 'Female', 'Male')))

df_y <- fread(Y_NCALLED_FILE, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE)

df <- merge(df, df_y, by='s')

colors <- pal_d3('category20')(20)[c(1,2)]
fills <- pal_d3('category20')(20)[c(11,12)]

create_pretty_cumulative(df, aes(impute_sex.f_stat), 'F-statistic', T_impute_sex,
                         xlim=c(-T_impute_sex,1.1), title='Cumulative Distribution of F-statistic', save_figure=TRUE, file=paste0(PLOTS,'05_F_stat_cdf'))

p <- ggplot(df, aes(x=impute_sex.f_stat, fill=imputed_sex)) +
  geom_histogram(binwidth=0.025, alpha=0.8, color='#7f7f7f') +
  scale_fill_manual(values=fills, limits=c('Male', 'Female')) +
  labs(x='X chromosome F-statistic',
       y='Count',
       title='',
       fill='Imputed Sex') +
  scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
  scale_y_continuous(label=scales::comma, breaks=scales::pretty_breaks(n=10)) +
  theme_minimal() +
  theme(axis.title.x=element_text(margin=ggplot2::margin(t=10)),
        axis.title.y=element_text(margin=ggplot2::margin(r=10)),
        plot.title=element_text(hjust=0.5)) +
  geom_vline(xintercept=T_impute_sex, linetype='dashed')

print(p)
ggsave(paste0(PLOTS, '05_imputesex_histogram', '.jpg'), p, width=160, height=90, units='mm')
ggsave(paste0(PLOTS, '05_imputesex_histogram', '.pdf'), p, width=160, height=90, units='mm')

#hashed out the below command as there is no missing sex info in data. change if not the case in other datasets
#if(any(df$phenotype.SEX == '' | df$phenotype.SEX == 'Not Reported' | df$phenotype.SEX == 'Unknown')) {
#  df$phenotype.SEX[df$phenotype.SEX == '' | df$phenotype.SEX == 'Not Reported' | df$phenotype.SEX == 'Unknown'] <- "Unknown"
#}

p <- ggplot(df, aes(x=impute_sex.f_stat, y=phenotype.SITE, colour=phenotype.SEX)) +
  geom_jitter(width=0, height=0.2, size=1, alpha=0.2, stroke=0.05) + 
  theme_minimal() +
  geom_vline(xintercept=T_impute_sex, linetype='dashed') +
  labs(x='X chromosome F-statistic',
       y='Location',
       color='Reported Sex') 

print(p)
ggsave(paste0(PLOTS, '05_imputesex_scatter_box', '.jpg'), p, width=160, height=90, units='mm')
ggsave(paste0(PLOTS, '05_imputesex_scatter_box', '.pdf'), p, width=160, height=90, units='mm')

df_false <- df %>% filter((impute_sex.f_stat > T_impute_sex & phenotype.SEX == 'F') | (impute_sex.f_stat < T_impute_sex & phenotype.SEX == 'M'))
df_false_plot <- df %>% filter(phenotype.SEX == 'Unknown' | (impute_sex.f_stat > T_impute_sex & phenotype.SEX == 'F') | (impute_sex.f_stat < T_impute_sex & phenotype.SEX == 'M'))

# Plots of gender estimates 
p <- ggplot(df, aes(x=impute_sex.f_stat, y=impute_sex.n_called, colour=phenotype.SEX)) +
  geom_point(size=0.5) + 
  labs(x='X chromosome F-statistic', y='Number of calls in Y', color='Reported Sex') +
  scale_color_d3('category10') +
  scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
  scale_y_continuous(label=scales::comma, breaks=scales::pretty_breaks(n=10)) +
  geom_point(data=df_false_plot, aes(x=impute_sex.f_stat, y=impute_sex.n_called), size=0.5) + 
  theme_minimal()
print(p)

ggsave(paste0(PLOTS, '05_imputesex_scatter', '.jpg'), p, width=160, height=90, units='mm')
ggsave(paste0(PLOTS, '05_imputesex_scatter', '.pdf'), p, width=160, height=90, units='mm')

df_out <- df_false %>% select(s)
write.table(df_out, file=SEXCHECK_LIST, quote=FALSE, row.names=FALSE, col.names=FALSE)
