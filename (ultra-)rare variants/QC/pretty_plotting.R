library(ggplot2)
library(dplyr)
library(ggsci)

create_pretty_forest <- function(df, title, save_figure=FALSE, file='file_out',
                                 width=160, height=90, scaling=1, y_label='', pass_order=NULL, print_p=TRUE,
                                 horizontal=TRUE, hline_at=0)
{
  if (!is.null(pass_order)) {
    df$label <- factor(df$label, levels = df$label[pass_order])
  }
  
  p <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
    geom_pointrange() + 
    geom_hline(yintercept=hline_at, lty=2)
  
  if (horizontal==TRUE) p <- p + coord_flip()
  
  p <- p + labs(title=title, y=y_label, x='') +
    theme_minimal()
  
  if (horizontal == TRUE) {
    p <- p + geom_text(aes(y=mean, label=p_vals), vjust=-1)
  } else {
    p <- p + geom_text(aes(y=mean, label=p_vals), hjust=-0.2)
  }
  
  if (print_p)  print(p)
  
  if (save_figure) {
    ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')
    ggsave(paste0(file, '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')
  }
  
  return(p)
}

create_pretty_hist <- function(df, aest, x_label, threshold=NULL, threshold_max=NULL,
                               file='file_out', title='', binwidth=0.002, width=160, height=90, scaling=1,
                               save_figure=FALSE, xlim=NULL, key_label='', print_p=TRUE)
{
  p <- ggplot(df, aest)
  geom_histogram(binwidth=binwidth, fill='#aec7e8', color='#1f77b4')
  ylim <- ggplot_build(p)$layout$panel_ranges[[1]]$y.range
  p <- ggplot(df, aest)
  
  if (is.null(aest$fill)) {
    p <- p + geom_histogram(binwidth=binwidth, fill='#aec7e8', color='#1f77b4')
  } else {
    p <- p + geom_histogram(binwidth=binwidth, color='grey50')
  }
  
  if (!is.null(xlim)) {
    p <- p + coord_cartesian(xlim=xlim)
  } 
  
  p <- p + labs(x=x_label, y='Count', title=title, fill=key_label) +
    scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
    scale_y_continuous(label=scales::comma, breaks=scales::pretty_breaks(n=10)) +
    theme_minimal() +
    theme(axis.title.x = element_text(margin=ggplot2::margin(t=10)),
          axis.title.y = element_text(margin=ggplot2::margin(r=10)),
          plot.title = element_text(hjust=0.5))
  
  if (!is.null(threshold)) {
    p <- p + geom_vline(xintercept=threshold, linetype='dashed') #+
    # annotate('text', label=paste0(x_label, ' = ', threshold), x=threshold, y=0.8*ylim[2], angle=90, vjust=-1)
  }
  if (!is.null(threshold_max)) {
    p <- p + geom_vline(xintercept=threshold_max, linetype='dashed') #+
    # annotate('text', label=paste0(x_label, ' = ', threshold_max), x=threshold_max, y=0.8*ylim[2], angle=90, vjust=2)
  }
  
  if (print_p) print(p)
  
  if (save_figure) {
    ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')
    ggsave(paste0(file, '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')
  }
  
  return(p)
}

create_pretty_density <- function(df, aest, x_label, threshold=NULL, threshold_max=NULL,
                                  file='file_out', title='', binwidth=0.002, width=160, height=90, scaling=1,
                                  save_figure=FALSE, xlim=NULL, key_label='', print_p=FALSE)
{
  p <- ggplot(df, aest)
  ylim <- ggplot_build(p)$layout$panel_ranges[[1]]$y.range
  p <- ggplot(df, aest)
  
  if (is.null(aest$fill)) {
    p <- p + geom_density(fill='#aec7e8', color='#1f77b4')
  } else {
    p <- p + geom_density(color='grey50')
  }
  
  if (!is.null(xlim)) {
    p <- p + coord_cartesian(xlim=xlim)
  } 
  
  p <- p + labs(x=x_label, y='Count', title=title, fill=key_label) +
    scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
    scale_y_continuous(label=scales::comma, breaks=scales::pretty_breaks(n=10)) +
    theme_minimal() +
    theme(axis.title.x = element_text(margin=ggplot2::margin(t=10)),
          axis.title.y = element_text(margin=ggplot2::margin(r=10)),
          plot.title = element_text(hjust=0.5))
  
  if (!is.null(threshold)) {
    p <- p + geom_vline(xintercept=threshold, linetype='dashed') #+
    # annotate('text', label=paste0(x_label, ' = ', threshold), x=threshold, y=0.8*ylim[2], angle=90, vjust=-1)
  }
  if (!is.null(threshold_max)) {
    p <- p + geom_vline(xintercept=threshold_max, linetype='dashed') #+
    # annotate('text', label=paste0(x_label, ' = ', threshold_max), x=threshold_max, y=0.8*ylim[2], angle=90, vjust=2)
  }
  
  if (print_p) print(p)
  
  if (save_figure) {
    ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')
    ggsave(paste0(file, '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')
  }
  
  return(p)
}

create_pretty_boxplots <- function(df, aes, aes_col, threshold=NULL,
                                   threshold_max=NULL, file='file_out', title='', x_label='', y_label='',
                                   key_label='', xlim=NULL, legend=FALSE, save_figure=FALSE, 
                                   width=160, height=90, scaling=1, facet=FALSE, facet_grid=NULL, jitter_size=0.5,
                                   outlier.shape=NA, n_ticks=10, print_p=FALSE, alpha=0.6)
{
  p = ggplot(df, aes) +
    geom_boxplot(outlier.shape=outlier.shape, coef=0, color='grey50', fill='grey95', show.legend=FALSE) + 
    geom_jitter(width=0.2, height=0, size=jitter_size, aes_col, show.legend=legend, alpha=alpha, stroke=0.05) + 
    coord_flip(ylim=xlim) +
    labs(title=title, x=y_label, y=x_label, color=key_label) + 
    scale_color_d3('category20') +
    scale_y_continuous(breaks = scales::pretty_breaks(n=n_ticks)) +
    guides(color = guide_legend(override.aes = list(size=2))) +
    theme_minimal() +
    theme(axis.title.x = element_text(margin = ggplot2::margin(t=10)),
          plot.title = element_text(hjust=0.5))
  
  if (!is.null(threshold)) {
    p <- p + geom_hline(yintercept=threshold, linetype='dashed')
  }
  if (!is.null(threshold_max)) {
    p <- p + geom_hline(yintercept=threshold_max, linetype='dashed')
  }
  if (facet){
    p <- p + facet_grid
  }
  
  if (print_p) print(p)
  
  if (save_figure) {
    ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')
    ggsave(paste0(file, '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')
  }
  return(p)
}

create_pretty_cumulative <- function(df, aes, x_label, threshold, threshold_max=NULL,
                                     file='file_out', title='', width=160, height=90, scaling=1, save_figure=FALSE,
                                     xlim=c(0,1), key_label='', print_p=FALSE)
{
  # These next ones are cdfs.
  p = ggplot(df, aes) + 
    stat_ecdf(geom='line', pad=FALSE) +
    geom_vline(xintercept=threshold, linetype='dashed') +
    coord_cartesian(xlim=xlim) +
    labs(x=x_label, y='Cumulative Proportion of Samples', title=title, color=key_label) +
    scale_color_d3('category10') +
    scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
    scale_y_continuous(breaks=scales::pretty_breaks(n=10)) +
    theme_minimal() +
    theme(axis.title.x = element_text(margin=ggplot2::margin(t=10)),
          axis.title.y = element_text(margin=ggplot2::margin(r=10)),
          plot.title = element_text(hjust=0.5))
  
  # p <- p + annotate('text', label=paste0(x_label, ' = ', threshold), x=threshold, y=mean(ylim), angle=90, vjust=-1)
  # ylim <- ggplot_build(p)$panel$ranges[[1]]$y.range
  # print(ylim)
  
  if (!is.null(threshold_max)) {
    p <- p + geom_vline(xintercept=threshold_max, linetype='dashed')# +
    # annotate('text', label=paste0(x_label, ' = ', threshold_max), x=threshold_max, y=mean(ylim), angle=90, vjust=2)
  }
  if (print_p) print(p)
  
  if (save_figure) {
    ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')
    ggsave(paste0(file, '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')
  }
  return(p)
}

create_pretty_scatter <- function(dt, aes, file='file_out', save_figure=FALSE, 
                                  key_label='', title='', limits=NULL, width=160, height=90, presentation=FALSE,
                                  add_final_layer=FALSE, final_layer=NULL, n_x_ticks=10, n_y_ticks=10, x_label=NULL,
                                  y_label=NULL, print_p=FALSE, gradient=FALSE, midpoint=0, gradient_title="")
{
  p <- ggplot(dt, aes) +
    geom_point(alpha=0.6, size=0.3) + 
    scale_color_d3('category20', limits=limits) +
    labs(title=title, color=paste0(key_label)) +
    scale_x_continuous(breaks=scales::pretty_breaks(n=n_x_ticks)) +
    scale_y_continuous(breaks=scales::pretty_breaks(n=n_y_ticks)) +
    theme_minimal() +
    theme(axis.title.x = element_text(margin=ggplot2::margin(t=10)),
          axis.title.y = element_text(margin=ggplot2::margin(r=10)),
          plot.title = element_text(hjust=0.5))
  
  if (add_final_layer) {
    cat("Adding final layer...\n")
    p <- p + guides(fill=FALSE) + geom_point(mapping=aes, data=final_layer) + scale_color_d3('category20')
  }
  
  if (gradient) {
    print(gradient_title)
    p <- p + scale_color_gradient2(
      low="blue", high="red", mid='grey50', midpoint=midpoint,
      name=gradient_title)
  }
  
  cat("Adding axis labels...\n")
  if (!is.null(x_label)) {
    p <- p + labs(x=x_label)
  }
  
  if (!is.null(y_label)) {
    p <- p + labs(y=y_label)
  }
  
  cat("Saving figure...\n")
  if (print_p) print(p)
  
  if (save_figure) {
    
    if (presentation == TRUE) {
      width <- 160
      height <- 90
      ggsave(paste0(file, '.jpg'), p, width=width, height=height, units='mm')
    } else {
      ggsave(paste0(file, '.jpg'), p, width=width, height=height, units='mm')
      ggsave(paste0(file, '.pdf'), p, width=width, height=height, units='mm')
    }
    
  }
  
  return(p)
}

create_pretty_qq_plot <- function(dt, aes, file='file_out', save_figure=FALSE,
                                  plot_title='', limits=NULL, width=110, height=110, n_x_ticks=10, n_y_ticks=10,
                                  x_label=TeX("$-\\log(p_{permutation})$"), 
                                  y_label=TeX("$-\\log(p_{observed})$"),
                                  n_to_include=NULL, cex_labels=1, print_p=TRUE, gradient=FALSE,
                                  gradient_title="", pval_col="pval")
{
  cat("Creating scatter-plot...\n")
  dt <- data.table(dt)
  setkeyv(dt, cols=pval_col)
  p <- create_pretty_scatter(dt, aes, file=file, save_figure=save_figure,
                             title=plot_title, limits=limits,
                             width=width, height=height, n_x_ticks=n_x_ticks, n_y_ticks=n_y_ticks,
                             x_label=x_label, y_label=y_label, gradient=gradient, gradient_title=gradient_title)
  
  cat("Adding y=x line...\n")
  p <- p + geom_abline(intercept=0, slope=1, color='indianred3') #+ coord_fixed()
  
  if (!is.null(n_to_include)) {
    cat("Adding labels...\n")
    p <- p + geom_label_repel(data=dt[(nrow(dt)-n_to_include+1):nrow(dt), ],
                              aes(label=labels), box.padding = 0.5, label.padding=0.1, point.padding = 0.2,
                              segment.color = 'grey50', size=cex_labels, segment.size=0.1)
  }
  
  cat("Created scatter-plot...\n")
  if (print_p) print(p)
  
  return(p)
}