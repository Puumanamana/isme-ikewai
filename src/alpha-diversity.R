source('loader.R')
source('util.R')

hue_cols <- c('PCA_Chem_Group', 'eruption')

##======================================##
##--- Load and preprocess metagenome ---##
##======================================##

metagenome <- load_metagenome()
metagenome <- filter_metagenome(metagenome, min_prevalence=2)

##======================================##
##----------- Compute metrics ----------##
##======================================##

metrics <- c("Observed", "InvSimpson", "Shannon", "Chao1")
data_all <- cbind(sample_data(metagenome),
                  estimate_richness(metagenome, measures=metrics))

# Add Pielou that's not implemented in Phyloseq
data_all$Pielou <- data_all$Shannon / log(data_all$Observed)
metrics <- c(metrics, 'Pielou')

##======================================##
##-------- Reshape data for plot -------##
##======================================##

data <- data_all[, c(hue_cols, metrics[metrics != 'Observed'])] %>%
  pivot_longer(names_to='metric', values_to='value', cols=-hue_cols) %>%
  pivot_longer(names_to='factor', values_to='level', cols=-c('metric', 'value'))

data$metric <- as.factor(data$metric)
data$factor <- as.factor(data$factor)
data$level <- factor(data$level, levels=levels(unlist(list(data_all[hue_cols]))))

label_levels <- sprintf(
  "%s\n(%s)", 
  rep(levels(data$metric), nlevels(data$factor)), 
  rep(levels(data$factor), each=nlevels(data$metric))
)

data$labels <- factor(sprintf("%s\n(%s)", data$metric, data$factor), levels=label_levels)

##======================================##
##------------ Plot and save -----------##
##======================================##

diversity_plot <- function(data=NULL, x=NULL, y=NULL, ncol=4) {
  p <- ggplot(data=data, aes_string(x=x[1], y=y, fill=x[length(x)])) +
    gginit() +
    geom_boxplot() +
    facet_wrap(~ labels, scales='free', ncol=ncol)
  return(p)
}

diversity_plot(data=data[data$factor=='eruption', ], x='level', y='value', ncol=4) +
  guides(fill = FALSE) + xlab('') + ylab('') +
  scale_fill_manual(values=colors$eruption)
ggsave(sprintf("%s/Figure-6 boxplot-alpha_diversity-eruption.pdf", io['figures']), width=12)

diversity_plot(data=data[data$factor=='PCA_Chem_Group', ], x='level', y='value', ncol=4) +
  guides(fill = FALSE) + xlab('') + ylab('') +
  scale_fill_manual(values=colors$PCA_Chem_Group)
ggsave(sprintf("%s/Figure-Sx boxplot-alpha_diversity-PCA_Chem_Group.pdf", io['figures']), width=12)

# diversity_plot(data, x=c('eruption', 'PCA_Chem_Group'), y=metrics) +
#   xlab('') +
#   theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
#   ggsave(sprintf("%s/boxplot-alpha_diversity-per-eruption-PCA_Chem_Group.pdf", io['figures']))
