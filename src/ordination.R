source('loader.R')
source('util.R')
library(vegan)

subsampling_level <- 5000

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
                            title="", size=2, cex=5) {
  
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
    p <- p + geom_point(aes_string(x=x, y=y, color=hue), shape=19, size=size, alpha=0.8)
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

##======================================##
##--- Load and preprocess metagenome ---##
##======================================##

metagenome <- load_metagenome()
metagenome <- filter_metagenome(metagenome, min_prevalence=2)

# Plot sample distribution for subsampling
plot_sample_distr(sample_sums(metagenome), cutoff=subsampling_level)
ggsave(sprintf("%s/sample_sizes.pdf", io['figures']), width=8, height=5)

# Subsampling
metagenome <- rarefy_even_depth(metagenome, subsampling_level, replace=FALSE, rngseed=42)

metadata <- as(sample_data(metagenome), 'data.frame')
hue_names <- c('PCA_Chem_Group', 'eruption')

##======================================##
##----------- CCA analysis -------------##
##======================================##

# Relative abundance
metagenome_relabund <- transform_sample_counts(metagenome, function(x) sqrt(x / sum(x)))

cca_model <- ordinate(metagenome_relabund, method='CCA', formula= ~ DO + pH + SPC + SO4 + NOx)

cca_components <- as.data.frame(scores(cca_model, display='sites'))
cca_components[, c(hue_names, 'site_id')] = metadata[, c(hue_names, 'site_id')]

cca_arrows <- as.data.frame(scores(cca_model, display="bp", scaling="species"))

for (hue_name in hue_names) {
  plot_components(data=cca_components, arrows=cca_arrows,
                  x='CCA1', y='CCA2', hue=hue_name,
                  label='site_id', size=7, cex=8) +
  scale_color_manual(values=colors[[hue_name]]) +
  
  if (hue_name == 'eruption') {
    ggsave(sprintf("%s/Figure-5 CCA-%s.pdf", io['figures'], hue_name), scale=2)
  } else {
    ggsave(sprintf("%s/Figure-S3 CCA-%s.pdf", io['figures'], hue_name), scale=2)
  }
}

##======================================##
##----------- NMDS analysis ------------##
##======================================##

method <- 'PCoA'

for (clade in names(clade_files)) {
  metagenome_clade <- filter_metagenome(metagenome, clade_file=clade_files[clade])
  metagenome_clade_relabund <- transform_sample_counts(metagenome_clade, function(x) sqrt(x / sum(x)))
  
  model <- ordinate(metagenome_clade_relabund, method=method, trymax=500, parallel=4)
  
  if (method == 'PCoA') {
    components <- as.data.frame(model$vectors[, 1:2])
  } else {
    components <- as.data.frame(scores(model))
  }
  colnames(components) <- sprintf('%s%s', method, 1:2)
  components[, hue_names] = metadata[rownames(components), hue_names]
  
  plot_components(data=components, x=paste0(method, '1'), y=paste0(method, '2'),
                  hue='PCA_Chem_Group', col='eruption') +
  scale_color_manual(values=colors$PCA_Chem_Group)
  
  for (hue_name in hue_names) {
    plot_components(data=components, # arrows=nmds_arrows, 
                    x=paste0(method, '1'), y=paste0(method, '2'), 
                    hue=hue_name) +
      
    scale_color_manual(values=colors[[hue_name]]) +
    ggsave(sprintf("%s/%s-%s-%s.pdf", io['figures'], method, hue_name, clade))
  }
}


