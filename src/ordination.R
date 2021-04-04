source('loader.R')
source('util.R')
library(vegan)
library(cowplot)

##======================================##
##---------- Plotting function ---------##
##======================================##

plot_sample_distr <- function(sizes, output=NULL, bins=30, cutoff=6000) {
  plot <- data.frame(frequency=sizes) %>% 
    ggplot(aes(frequency)) +
    gginit() +
    geom_histogram(fill='#006daa', color='black', bins=bins) +
    scale_x_log10() +
    xlab('Sample size') + ylab('Frequency') +
    geom_vline(xintercept=c(cutoff), color="red", linetype=2) +
    annotate("text", x=cutoff, y=7, label=cutoff, color='red', fontface='bold')
  
  return(plot)
}

plot_components <- function(data=NULL, x=NULL, y=NULL, hue=NULL, arrows=NULL, col=NULL, label=NULL,
                            title="", size=2, cex=5, alpha=0.8) {
  
  p <- ggplot(data) +
    gginit() + 
    coord_fixed() +
    geom_vline(xintercept=c(0), color="grey70", linetype=2) +
    geom_hline(yintercept=c(0), color="grey70", linetype=2) + 
    labs(title=title)
  
  # Plot textual labels instead of circles
  if (!is.null(label)) {
    p <- p + geom_text(aes_string(x=x, y=y, color=hue, label=label), size=size, fontface = "bold")
  } else {
    p <- p + geom_point(aes_string(x=x, y=y, color=hue), shape=19, size=size, alpha=alpha)
  }
  
  # Overlay arrows
  if (!is.null(arrows)) {
    arrows_scale <- 0.4 * max(data[, c(x,y)]) / max(arrows)
    arrows_rescale <- arrows * arrows_scale
    label_nudge <- max(arrows_rescale) * 0.2
    
    p <- p + 
      geom_segment(data=arrows_rescale, aes_string(x=0, y=0, xend=x, yend=y), 
                   arrow=arrow(length=unit(0.2, "cm"))) + 
      geom_text(data=arrows_rescale, aes_string(x=x, y=y), label=rownames(arrows_rescale),
                cex=cex, nudge_x=label_nudge, nudge_y=label_nudge, fontface='bold')
  } else if (!is.null(col)) {
    # Can't do both arrows and facets. If no arrows, facet_grid is possible
    p <- p + facet_wrap(as.formula(sprintf('~ %s', col)))
  }
  
  return(p)
  
}

##==================##
##--- Parameters ---##
##==================##

mode <- 'single'
subsampling_level <- c(merged=8000, single=3000)
seed <- 123
set.seed(seed)

##======================================##
##--- Load and preprocess metagenome ---##
##======================================##

metagenome <- load_metagenome(abundance_file=sprintf('../data/%s/abundance.csv',  mode),
                              metadata_file=sprintf('../data/%s/metadata.csv', mode))
# metadata <- sample_data(metagenome)
# sample_data(metagenome) <- metadata
metagenome <- subset_samples(metagenome, Site_id %in% 5:20)

# Plot sample distribution for subsampling
plot_sample_distr(sample_sums(metagenome), cutoff=subsampling_level[mode])
ggsave(sprintf("%s/sample_sizes.pdf", io['figures']), width=8, height=5)

metagenome_filt <- filter_prevalence(metagenome, min_prevalence=2)
metagenome_subsampled <- rarefy_even_depth(metagenome_filt, subsampling_level[mode], replace=FALSE, rngseed=seed)
metagenome_relabund <- transform_sample_counts(metagenome_subsampled, function(x) sqrt(x / sum(x)))

# Subsampling
metadata <- as(sample_data(metagenome), 'data.frame')

##======================================##
##----------- CCA analysis -------------##
##======================================##

cca_model <- ordinate(metagenome_relabund, method='CCA', formula= ~ DO + pH + SPC + SO4 + NOx)

cca_components <- as.data.frame(scores(cca_model, display='sites'))
cca_components[, colnames(metadata)] = metadata[rownames(cca_components),]

cca_arrows <- as.data.frame(scores(cca_model, display="bp", scaling="species"))

plot_components(data=cca_components, arrows=cca_arrows, x='CCA1', y='CCA2', 
                hue='Eruption', label='Site_id', size=5, cex=8) +
    scale_color_manual(values=colors[['Eruption']])
ggsave(sprintf("%s/Figure-5 CCA-eruption.pdf", io['figures']), scale=2)

plot_components(data=cca_components, arrows=cca_arrows, x='CCA1', y='CCA2', 
                hue='group', label='Site_id', size=7, cex=8)
ggsave(sprintf("%s/Figure-S3 CCA-groups.pdf", io['figures']), scale=2)

# plot_components(data=cca_components, arrows=cca_arrows, x='CCA1', y='CCA2', 
#                 hue='PCA_Grp', label='Site_id', size=7, cex=8) +
#   scale_color_manual(values=colors[['PCA_Grp']])
# ggsave(sprintf("%s/Figure-S3 CCA-PCA_Grp.pdf", io['figures']), scale=2)

##======================================##
##----------- NMDS analysis ------------##
##======================================##

run_ordination <- function(mg, method='NMDS', dist='bray', trymax=500) {
  
  if (method == 'PCoA') {
    model <- ordinate(mg, method='PCoA', distance=dist)
    components <- as.data.frame(model$vectors[, 1:2])
  } else {
    if (endsWith(dist, 'unifrac')) {
      distances <- UniFrac(mg, weighted=dist[1]=='w', normalized=TRUE, 
                           parallel=FALSE, fast=TRUE)
      model <- metaMDS(distances, trymax=trymax)
    } else {
      model <- ordinate(mg, method='NMDS', distance=dist, trymax=trymax)
    }
    components <- as.data.frame(model$points)
  }
  
  meta <- sample_data(mg)[rownames(components),]
  meta[, sprintf('PC%s', 1:2)] <- components
  
  return(meta)
}

method <- 'NMDS'
dist <- 'bray'

components <- run_ordination(metagenome_subsampled, method=method, dist=dist, trymax=500)

# NMDS colored by eruption
plot_components(data=components, alpha=1, size=3, x='PC1', y='PC2', hue='Eruption') + 
  stat_ellipse(geom="polygon", aes(x=PC1, y=PC2, fill=Eruption), alpha=0.2, show.legend=FALSE, level=0.9)
ggsave(sprintf("%s/%s-eruption.pdf", io['figures'], method))

# NMDS colored by Year-Month
plot_components(data=components, alpha=1, size=3, x='PC1', y='PC2', hue='YM') + 
  stat_ellipse(geom="polygon", aes(x=PC1, y=PC2, fill=YM), alpha=0.2, show.legend=FALSE, level=0.5)
ggsave(sprintf("%s/%s-YM.pdf", io['figures'], method))

# NMDS colored by location
plot_components(data=components, alpha=1, size=5, x='PC1', y='PC2', hue='group', label='Site_id')

plots <- list()
subsampling <- data.frame(
  single=c('2017-11'=3000, '2018-03'=10000, '2018-08'=5000, '2018-11'=4000, '2019-03'=7000),
           # 'PreEruption'=6000, 'PostEruption'=2000),
  merged=c('2017-11'=3000, '2018-03'=15000, '2018-08'=10000, '2018-11'=10000, '2019-03'=12000)
           # 'PreEruption'=12000, 'PostEruption'=12000)
)

for (ym in rownames(subsampling)) {
  # ym <- rownames(subsampling)[5]
  mg_sub <- subset_samples(metagenome, YM==ym)
  # mg_sub <- subset_samples(mg_sub, group != '16-20')
  summary <- data.frame(
    site=metadata[sample_names(mg_sub),'Site_id'], 
    size=sample_sums(mg_sub)
  )
  summary <- summary[order(summary$size),]
  plot_sample_distr(summary$size, cutoff=subsampling[ym, mode], bins=20)
  
  mg_sub <- rarefy_even_depth(mg_sub, subsampling[ym, mode], replace=FALSE, rngseed=seed)
  # mg_sub <- transform_sample_counts(mg_sub, function(x) asin(sqrt(x / sum(x))))
  components <- run_ordination(mg_sub, method=method, dist=dist)
  plots[[ym]] <- plot_components(data=components, x='PC1', y='PC2', hue='group', label='Site_id', size=5)
}

combined_plot <- plot_grid(
    plots[['2017-11']] + theme(legend.position="none") + labs(title='2017-11'), 
    plots[['2018-03']] + theme(legend.position="none") + labs(title='2018-03'),
    plots[['2018-08']] + theme(legend.position="none") + labs(title='2018-08'),
    plots[['2018-11']] + theme(legend.position="none") + labs(title='2018-11'),
    plots[['2019-03']] + labs(title='2019-03'),
    ncol=3, align = "hv", rel_widths = c(8,8)
)
combined_plot
ggsave(sprintf("%s/NMDS-group-by-YM.pdf", io['figures']), combined_plot, height=11, width=17)
  