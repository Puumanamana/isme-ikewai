source('loader.R')
source('plot_utils.R')
library(vegan)

fig_dir <- "../figures"

##======================================##
##--- Load and preprocess metagenome ---##
##======================================##

metagenome <- load_metagenome()
metagenome <- filter_metagenome(metagenome, min_prevalence=1)
metagenome <- transform_sample_counts(metagenome, function(x) sqrt(x / sum(x)))

metadata <- as(sample_data(metagenome), 'data.frame')
hue_names <- c('PCA_Chem_Group', 'eruption')

##======================================##
##----------- CCA analysis -------------##
##======================================##

cca_model <- ordinate(metagenome, method='CCA', 
                      formula= ~ DO + pH + SPC + SO4 + NOx)

cca_components <- as.data.frame(scores(cca_model, display='sites'))
cca_components[, hue_names] = metadata[, hue_names]

cca_arrows <- as.data.frame(scores(cca_model, display="bp", scaling="species"))

for (hue_name in hue_names) {
  plot_components_with_arrows(data=cca_components, arrows=cca_arrows,
                              x='CCA1', y='CCA2', hue=hue_name,
                              title='Canonical Correspondance Analysis')
  ggsave(sprintf("%s/CCA-%s.pdf", fig_dir, hue_name))
}

##======================================##
##----------- NMDS analysis ------------##
##======================================##

nmds_model <- metaMDS(otu_table(metagenome), trymax=500, parallel=4)
# nmds_model <- ordinate(metagenome_normalized, method='NMDS', trymax=500, parallel=4)
nmds_env <- envfit(nmds_model ~ DO + pH + SPC + SO4 + NOx, data=metadata, perm=999)

nmds_components <- as.data.frame(scores(nmds_model))
nmds_components[, hue_names] = metadata[, hue_names]

nmds_arrows <- as.data.frame(nmds_env$vectors$arrows)

for (hue_name in hue_names) {
  plot_components_with_arrows(data=nmds_components, arrows=nmds_arrows, 
                              x='NMDS1', y='NMDS2', hue=hue_name,
                              title='Non-metric Multidimentional Scaling')
  ggsave(sprintf("%s/NMDS-%s.pdf", fig_dir, hue_name))
}
