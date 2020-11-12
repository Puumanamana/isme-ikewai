#============ Set up environment ============#

source('loader.R')
source('util.R')

library(DESeq2)
library(vegan)
library(BiocParallel)

pval <- 0.05
hue_cols <- c('PCA_Chem_Group', 'eruption')

threads <- 3
register(MulticoreParam(threads))

#=========== Loading and Filtering ===========#

metagenome <- load_metagenome()
metagenome <- subset_samples(metagenome, PCA_Chem_Group %in% c('C', 'D'))
metagenome <- speedyseq::tax_glom(metagenome, 'Genus')
metagenome <- filter_metagenome(metagenome, min_prevalence=1)

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

#========== Filter and sort result ============#

sig_otu_table <- results(diagdds, alpha=pval, pAdjustMethod='BH',parallel=TRUE)
sig_otu_table <- sig_otu_table[order(sig_otu_table$padj, na.last=NA), ]
sig_otu_table$padj <- signif(sig_otu_table$padj, 3)

sig_otu_table <- cbind(
  as(sig_otu_table, "data.frame"), 
  as(tax_table(metagenome)[rownames(sig_otu_table), ], "matrix")
)

#==================================================#
#------------ Plot OTU distributions --------------#
#==================================================#

taxa_to_plot <- list(dsrB=c('Desulfurivibrio', 'Desulfovibrio', 'Aeromonas', 'Alcanivorax', 
                            'Hymenobacter', 'Desulfobulbus'), 
                     nirS=c('Pseudomonas', 'Aquabacterium', 'Hydrogenophaga', 'Halomonas', 
                            'Desulfurivibrio', 'Bradyrhizobium', 'Pseudoxanthomonas'))

data <- c()
for (clade in names(taxa_to_plot)) {
  mg_subset <- subset_taxa(metagenome, Genus %in% taxa_to_plot[[clade]])
  abundance <- as.data.frame(otu_table(mg_subset))
  colnames(abundance) <- tax_table(mg_subset)[colnames(abundance), 'Genus']

  data <- abundance %>%
    bind_cols(as(sample_data(mg_subset)[, hue_cols], 'data.frame')) %>%
    pivot_longer(-hue_cols, names_to='Genus', values_to='abundance') %>%
    pivot_longer(-c('Genus', 'abundance'), names_to='factor', values_to='level') %>%
    left_join(sig_otu_table[, c('Genus', 'padj')]) %>%
    mutate(clade=clade, label=sprintf('%s (%s, padj=%s)', Genus, clade, padj)) %>%
    bind_rows(data) 
}

ggplot(data[data$factor=='eruption', ], aes(x=Genus, y=abundance, fill=level)) +
  gginit() + 
  geom_boxplot() +
  scale_fill_manual(values=colors$eruption) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base=10)) +
  facet_wrap(' ~ label', scale="free") +
  labs(x='', y='abundance')
ggsave(sprintf('%s/Figure-8A distr_eruption.pdf', io['figures']), width=50, height=30, units="cm")

ggplot(data[data$factor=='PCA_Chem_Group', ], aes(x=Genus, y=abundance, fill=level)) +
  gginit() + 
  geom_boxplot() +
  scale_fill_manual(values=colors$PCA_Chem_Group) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base=10)) +
  facet_wrap(' ~ label', scale="free") +
  labs(x='', y='abundance')
ggsave(sprintf('%s/Figure-8B distr_PCA_Chem_Group.pdf', io['figures']), width=50, height=30, units="cm")

