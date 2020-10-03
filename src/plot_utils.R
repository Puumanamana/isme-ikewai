library(ggplot2)
library(ggsci)
library(dplyr)

gginit <- function() {
  return(
    theme_linedraw() + 
    theme(
        title=element_text(size=16),
        strip.text.x=element_text(size=14),
        panel.grid=element_line("white"),
        axis.line=element_line("gray25"),
        axis.text=element_text(size=12, color="gray25"),
        axis.title=element_text(size=12, color="gray25"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12)
    )
  )
}

plot_components_with_arrows <- function(data=NULL, x=NULL, y=NULL, hue=NULL, arrows=NULL, 
                                        title="", size=2, cex=4) {
  pal <- pal_jco()(nlevels(data[[hue]]))
  
  arrows_scale <- 0.4 * max(data[, c(x,y)]) / max(arrows)
  arrows_rescale <- arrows * arrows_scale
  label_nudge <- max(arrows_rescale) * 0.2
  
  p <- ggplot(data) +
    gginit() + 
    geom_point(aes_string(x=x, y=y, color=hue), shape=19, size=size, alpha=0.8) +
    coord_fixed() +
    scale_color_manual(values=pal) +
    geom_segment(data=arrows_rescale, aes_string(x=0, y=0, xend=x, yend=y), 
                 arrow=arrow(length=unit(0.2, "cm"))) +
    geom_vline(xintercept=c(0), color="grey70", linetype=2) +
    geom_hline(yintercept=c(0), color="grey70", linetype=2) +
    geom_text(data=arrows_rescale, aes_string(x=x, y=y), label=rownames(arrows_rescale),
              cex=cex, nudge_x=label_nudge, nudge_y=label_nudge) +
    labs(title=title)
  
  return(p)
  
}

stacked_barplot <- function(metagenome, x, hue='Phylum') {

  mg <- tax_glom(metagenome, hue) # Group OTUs by `hue`
  
  # Group by multiple columns
  factor_data <- sapply(x, function(y) get_variable(mg, y))
  sample_data(mg)$combined <- apply(factor_data, 1, paste0, collapse=';')
  mg <- merge_samples(mg, 'combined')
  
  multiindex <- t(as.data.frame(strsplit(sample_names(mg), ';')))
  # Fix merge_samples metadata (GitHub issue #243)
  for (i in 1:length(x)) {
    sample_data(mg)[[x[i]]] <- factor(
      as.character(multiindex[, i]), levels=levels(sample_data(metagenome)[[x[i]]])
    )
  }
  
  mg <- transform_sample_counts(mg, function(y) y/sum(y))
  data <- psmelt(mg)

  rank_sums <- data %>% 
    group_by_at(hue) %>% 
    summarize(Abundance=sum(Abundance)) %>% 
    arrange(Abundance)
  
  data[[hue]] <- factor(data[[hue]], levels=rank_sums[[hue]])
  data <- data[order(data[[hue]]), ]
  
  nhue <- nlevels(data[[hue]])
  if (nhue < 52) {
    palette <- pal_igv()(nhue)
  } else {
    palette <- rep(pal_igv()(51), 1+int(nhue/51))[1:nhue]
  }
  
  p <- ggplot(data=data, aes_string(x=x[1], y='Abundance', fill=hue)) +
    gginit() + 
    geom_bar(aes(), stat="identity", position="stack", color='black', 
             size=0.5, width=0.7) +
    scale_fill_manual(values=palette) +
    ylab('Proportion') + 
    labs(title=sprintf('Taxonomic composition (%s)', hue))

  if (length(x) > 1) {
    p <- p + facet_grid(as.formula(sprintf('~ %s', x[-1])))
  }
  
  return(p)
}


diversity_plot <- function(data=NULL, x=NULL, y=NULL) {
  
  df <- reshape2::melt(data[c(x, y)], id_vars=x, variable.name='metric', value.name='score')
  p <- ggplot(data=df, aes_string(x=x[1], y='score', fill=x[length(x)])) +
    gginit() + 
    geom_boxplot() +
    facet_wrap( . ~ metric, scales='free_y')
  return(p)
}