source('loader.R')
source('plot_utils.R')

fig_dir <- "../figures"

##======================================##
##--- Load and preprocess metagenome ---##
##======================================##

metagenome <- load_metagenome()
metagenome <- filter_metagenome(metagenome, min_prevalence=1)

##======================================##
##----------- Compute metrics ----------##
##======================================##

metrics <- c("Observed", "InvSimpson", "Shannon", "Chao1")
data <- cbind(sample_data(metagenome),
              estimate_richness(metagenome, measures=metrics))

##======================================##
##------------ Plot and save -----------##
##======================================##

diversity_plot(data, x=c('eruption'), y=metrics) +
  ggsave(sprintf("%s/boxplot-alpha_diversity-per-eruption.pdf", fig_dir))
diversity_plot(data, x=c('PCA_Chem_Group'), y=metrics) +
  ggsave(sprintf("%s/boxplot-alpha_diversity-per-PCA_Chem_Group.pdf", fig_dir))
diversity_plot(data, x=c('eruption', 'PCA_Chem_Group'), y=metrics) +
  ggsave(sprintf("%s/boxplot-alpha_diversity-per-eruption-PCA_Chem_Group.pdf", fig_dir))
