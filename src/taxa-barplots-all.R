source('util.R')
source('loader.R')
source('~/code/r-util/taxonomy_plot.R')
library(cowplot)

mode <- 'single'
outdir <- '~/code/isme-ikewai/figures'
  
##======================================##
##---------- Load metagenome -----------##
##======================================##

metagenome <- load_metagenome(abundance_file=sprintf('../data/%s/abundance.csv',  mode),
                              metadata_file=sprintf('../data/%s/metadata.csv', mode))
metagenome <- subset_samples(metagenome, Site_id %in% 5:20)

##=====================================##
##---------- Barplots (all) -----------##
##=====================================##

taxa_barplot(metagenome, 'Phylum', 'group', thresh=1e-2) +
  guides(fill=guide_legend(nrow=30, reverse=T))
ggsave(sprintf('%s/barplot-phylum-site.pdf', outdir), height=6, width=6)

taxa_barplot(metagenome, 'Class', 'YM', thresh=1e-2) +
  guides(fill=guide_legend(nrow=30, reverse=T))
ggsave(sprintf('%s/barplot-phylum-YM.pdf', outdir), height=6, width=6)

taxa_barplot(metagenome, 'Phylum', c('Site_id', 'Year', 'Month'), thresh=1e-2) +
  guides(fill=guide_legend(nrow=35, reverse=T)) +
  theme(panel.grid.major=element_line(size=0.05), panel.grid.minor=element_blank())
ggsave(sprintf('%s/barplot-phylum-site-year-month.pdf', outdir), height=6, width=6)

mg_proteo <- subset_taxa(metagenome, Phylum == 'Proteobacteria')
taxa_barplot(mg_proteo, 'Order', c('Site_id', 'Year', 'Month'), thresh=1e-2) +
  guides(fill=guide_legend(nrow=35, reverse=T))
ggsave(sprintf('%s/barplot-proteo_order-site-year-month.pdf', outdir), height=6, width=6)

##==========================================##
##---------- Barplots (dsr+nirS) -----------##
##==========================================##

mg_dsr <- subset_metagenome(metagenome, io[['dsrB']])
dsr_plot <- taxa_barplot(mg_dsr, 'Genus', 'group', thresh=1e-3) +
  guides(fill=guide_legend(nrow=20, reverse=T))

mg_nirS <- subset_metagenome(metagenome, io[['nirS']])
nirs_plot <- taxa_barplot(mg_nirS, 'Genus', 'group', thresh=1e-3) +
  guides(fill=guide_legend(nrow=20, reverse=T))

combined_plot <- plot_grid(
  dsr_plot + labs(title='Sulfate-reducing bacteria'),
  nirs_plot + labs(title='Denitrifying bacteria'),
  ncol=1, align = "hv"
)
combined_plot
ggsave(sprintf("%s/barplot-nirS-dsrB_by-genus.pdf", outdir), combined_plot, height=11, width=6)

