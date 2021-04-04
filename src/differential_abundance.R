#============ Set up environment ============#

source('loader.R')
source('util.R')

library(DESeq2)
library(vegan)
library(BiocParallel)

pval <- 0.05
hue_cols <- c('group', 'YM')
mode <- 'single'

threads <- 3
register(MulticoreParam(threads))

#=========== Loading and Filtering ===========#

metagenome <- load_metagenome(abundance_file=sprintf('../data/%s/abundance.csv',  mode),
                              metadata_file=sprintf('../data/%s/metadata.csv', mode))
metagenome <- subset_samples(metagenome, Site_id %in% 5:20)
# metagenome <- subset_samples(metagenome, PCA_Grp %in% c('C', 'D'))
metagenome <- speedyseq::tax_glom(metagenome, 'Genus')
metagenome <- filter_prevalence(metagenome, min_prevalence=2)

#==================================================#
#---------------- DESeq 2 package -----------------#
#==================================================#

## Script from tutorial at http://joey711.github.io/phyloseq-extensions/DESeq2.html

#============= Fit DESeq2 model ===============#
design <- as.formula(sprintf(" ~ %s", paste(hue_cols, collapse='+')))

diagdds <- phyloseq_to_deseq2(metagenome, design)

# calculate geometric means prior to estimate size factors
# DESeq2 cannot normalize when every OTU contains at least one zero
geom_mean <- function(x) { exp(mean(log(x[x>0]), na.rm=TRUE)) }

geoMeans <- apply(counts(diagdds), 1, geom_mean)
diagdds <- estimateSizeFactors(diagdds, geoMeans=geoMeans)
diagdds <- DESeq(diagdds, fitType="local", test="Wald", parallel=TRUE, quiet=TRUE)
# diagdds <- DESeq(diagdds, parallel=TRUE, quiet=TRUE)

#========== Filter and sort result ============#
sig_otu_table <- results(diagdds, alpha=pval, pAdjustMethod='BH',parallel=TRUE)
sig_otu_table <- sig_otu_table[order(sig_otu_table$padj, na.last=NA), ]
sig_otu_table$padj <- signif(sig_otu_table$padj, 3)
sig_otu_table <- sig_otu_table[sig_otu_table$padj<0.5,]

sig_otu_table <- cbind(
  as(sig_otu_table, "data.frame"), 
  as(tax_table(metagenome)[rownames(sig_otu_table), ], "matrix")
)

sig_otu_table[, c('log2FoldChange', 'padj', 'Family', 'Genus')] %>% 
  head(50)

#==================================================#
#------------ Plot OTU distributions --------------#
#==================================================#

# filtering function - turns outliers into NAs to be removed
filter_lims <- function(x){
  l <- boxplot.stats(x)$stats[1]
  u <- boxplot.stats(x)$stats[5]
  
  for (i in 1:length(x)){
    x[i] <- ifelse(x[i]>l & x[i]<u, x[i], NA)
  }
  return(x)
}

taxa_to_plot <- list(dsrB=c('Desulfurivibrio', 'Desulfovibrio', 'Aeromonas', 'Alcanivorax', 
                            'Hymenobacter', 'Desulfobulbus', 'Candidatus_Desulforudis',
                            'Chromatiaceae_unclassified', 'Desulfobulbaceae', 'Thiobacillus'), 
                     nirS=c('Pseudomonas', 'Aquabacterium', 'Hydrogenophaga', 'Halomonas', 
                            'Desulfurivibrio', 'Bradyrhizobium', 'Pseudoxanthomonas', 'Denitratisoma'),
                     methane=c('Methanobacterium', 'Methyloparacoccus', 'Methanosaeta',
                               'Methylophilaceae_unclassified', 'Methylotenera', 'Methanoregula', 
                               'UBA6140'))
data <- c()
for (clade in names(taxa_to_plot)) {
  mg_subset <- subset_taxa(metagenome, 
                           (Genus %in% taxa_to_plot[[clade]]) | (Family %in% taxa_to_plot[[clade]]))

  data <- psmelt(mg_subset) %>%
    inner_join(sig_otu_table[, c('Genus', 'padj')]) %>%
    mutate(clade=clade, label=sprintf('%s (%s, padj=%s)', Genus, clade, padj)) %>%
    bind_rows(data)
}

data %>%
  # group_by(label, group, YM) %>%
  # mutate(Abundance=filter_lims(Abundance)) %>%
  ggplot(aes(x=Eruption, y=Abundance, fill=group)) +
  gginit() + 
  geom_boxplot() +
  scale_fill_manual(values=colors$group) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base=10)) +
  facet_wrap(' ~ label', scale="free") +
  labs(x='', y='Abundance')
ggsave(sprintf('%s/Figure-8B distr_PCA_Grp.pdf', io['figures']), width=50, height=30, units="cm")

data %>%
  # group_by(label, Eruption) %>%
  # mutate(Abundance=filter_lims(Abundance)) %>%
  ggplot(aes(x=Genus, y=Abundance, fill=Eruption)) +
  gginit() +
  geom_boxplot(na.rm = TRUE, coef=5) +
  scale_fill_manual(values=colors$Eruption) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base=10)) +
  facet_wrap('~ label', scale="free") +
  labs(x='', y='Abundance')
ggsave(sprintf('%s/Figure-8A distr_eruption.pdf', io['figures']), width=50, height=30, units="cm")

