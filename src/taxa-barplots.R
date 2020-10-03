source('loader.R')
source('plot_utils.R')

fig_dir <- "../figures"

##======================================##
##---------- Load metagenome -----------##
##======================================##

metagenome <- load_metagenome()

##======================================##
##------- Custom stacked barplot -------##
##======================================##

for (clade in c('nirS', 'dsrB')) {
  mg <- filter_metagenome(
    metagenome, min_prevalence=0, 
    clade_file=sprintf('../data/functional-clades/%s_genuses.tsv', clade)
  )
  
  for (rank in c('Genus')) {
    stacked_barplot(mg, 'eruption', rank)
    ggsave(sprintf("%s/%s-barplot-eruption-per-%s.pdf", fig_dir, clade, rank))
    
    stacked_barplot(mg, 'PCA_Chem_Group', rank) 
    ggsave(sprintf("%s/%s-barplot-PCA_Chem_Group-per-%s.pdf", fig_dir, clade, rank))
  }
}

##======================================##
##---------- Other cool stuff ----------##
##======================================##

# plot_heatmap(metagenome, sample.label='site_id', sample.order='site_id', taxa.label="Genus", method="NMDS", distance="bray")

# plot_net(metagenome, distance="(A+B-2*J)/(A+B)", type="samples",color="PCA_Chem_Group", maxdist=0.7)
# plot_net(metagenome, distance="(A+B-2*J)/(A+B)", type="taxa", color="Phylum", maxdist=0.7)
