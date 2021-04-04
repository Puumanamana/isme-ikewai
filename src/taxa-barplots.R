source('loader.R')
source('util.R')
library(cowplot)

factors <- c('Site_id', 'Eruption')
multilevel <- TRUE
mode <- 'single'
rank <- 'Genus'

##======================================##
##---------- Load metagenome -----------##
##======================================##

metagenome <- load_metagenome(abundance_file=sprintf('../data/%s/abundance.csv',  mode),
                              metadata_file=sprintf('../data/%s/metadata.csv', mode))
metagenome <- subset_samples(metagenome, Site_id %in% 5:20)
# metagenome <- subset_samples(metagenome, PCA_Grp %in% c('C', 'D'))
metagenome <- speedyseq::tax_glom(metagenome, rank)

# Combine hue_cols for multilevel plot
sample_data(metagenome)$combined <- factor(sprintf(
  "%s;%s", get_variable(metagenome, factors[1]), get_variable(metagenome, factors[2])
))

##======================================##
##-------------- Format data -----------##
##======================================##

merge_and_format_samples <- function(mg, f) {
  # Need for a specific function to handle (GitHub issue #243)
  mg_ <- merge_samples(mg, f)
  sample_data(mg_)[[f]] <- factor(sample_names(mg_))
  mg_ <- transform_sample_counts(mg_, function(y) y/sum(y))
  
  # Long format
  data <- psmelt(mg_)[, c(rank, 'Abundance', f)] %>%
    rename(level=f) %>%
    mutate(factor=all_of(f))
  
  return(data)
}

data_all <- list()

for (clade in names(io)[2:3]) {
  if (is.na(clade)) {
    clade <- 'all'
    metagenome_clade <- metagenome
  } else {
    metagenome_clade <- filter_metagenome(metagenome, clade_file=io[clade])
  }

  # Aggregate, normalize and melt data for each factor
  hues <- factors
  if (multilevel) {
    hues <- 'combined'
  }
  
  data <- c()
  for (hue in hues) {
    tmp <- merge_and_format_samples(metagenome_clade, hue)
    tmp$level <- factor(tmp$level, levels=levels(get_variable(metagenome, hue)))
    data <- data %>% bind_rows(tmp)
  }
  
  if(multilevel) {
    data <- data %>% 
      separate(level, into=factors, sep=';')
    for (f in factors) {
      data[[f]] <- factor(data[[f]], levels=levels(get_variable(metagenome, f)))
    }
  }
  
  rank_sums <- data %>% 
    group_by_at(rank) %>% 
    summarize(Abundance=sum(Abundance)) %>% 
    arrange(Abundance)
  
  data[[rank]] <- factor(data[[rank]], levels=rank_sums[[rank]])
  data_all[[clade]] <- data[order(data[[rank]]), ]
}


##===================================##
##------ Custom stacked barplot -----##
##===================================##

plots <- list()
for (clade in names(io)[2:3]) {
  data <- data_all[[clade]]
  data <- data[data$Abundance > 0,]
  
  ntaxa <- nlevels(data[[rank]])
  if (ntaxa < 52) {
    palette <- pal_igv()(ntaxa)
  } else {
    palette <- rep(pal_igv()(51), 1+floor(ntaxa/51))[1:ntaxa]
  }

  if (multilevel) {
    plots[[clade]] <- 
      ggplot(data=data, aes_string(x=factors[1], y='Abundance', fill=rank)) +
      gginit() +
      geom_bar(stat="identity", position="stack", color='black', size=0.5, width=0.7) +
      scale_fill_manual(values=palette) +
      ylab('Proportion') + xlab('') +
      guides(fill=guide_legend(title=sprintf("%s for %s", rank, clade), reverse=TRUE)) + 
      facet_wrap(Eruption ~ .)
      
    
  } else {
    plots[[clade]] <- 
      ggplot(data=data, aes_string(x='level', y='Abundance', fill=rank)) +
      gginit() + 
      geom_bar(stat="identity", position="stack", color='black', 
               size=0.5, width=0.7) +
      scale_fill_manual(values=palette) +
      ylab('Proportion') + xlab('') +
      guides(fill=guide_legend(title=sprintf("%s for %s", rank, clade), reverse=TRUE)) + 
      facet_grid(~ factor, scales='free')
  }

}

plot <- plot_grid(
  plot_grid(
    plots[['nirS']] + theme(legend.position = "none"), 
    plots[['dsrB']] + theme(legend.position = "none"),
    ncol=1, align = "hv"
  ), plot_grid(
    get_legend(plots[['nirS']]),
    get_legend(plots[['dsrB']]), ncol=1)
  , rel_widths = c(8,3)
)

if(multilevel){
  ggsave(sprintf("%s/Figure-S4 taxa-barplot.pdf", io['figures']), plot, height=11, width=17)
} else{
  ggsave(sprintf("%s/Figure-7 taxa-barplot.pdf", io['figures']), plot, height=11, width=17)
}

