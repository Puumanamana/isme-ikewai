##==================================================##
##----------- Metagenome loading utils -------------##
##==================================================##

library(phyloseq)
library(data.table)

##============ txt file loader into phyloseq object ============##

fast_table_load <- function(filename, threads=1, drop=c(), row.names=1) {
  data <- fread(filename, nThread=threads, drop=drop, header=T, blank.lines.skip=T)
  data <- as.data.frame(data)
  
  if (row.names > 0) {
    rownames(data) <- data[[row.names]]
    data[, row.names] <- NULL
  }
  return(data)
}

load_metagenome <- function(abundance_file='../data/abundance.csv', 
                            metadata_file='../data/metadata.csv',
                            taxonomy_file='../data/taxonomy.csv',
                            tree_file='../data/tree.nwk',
                            threads=4) {
  abundance <- fast_table_load(abundance_file, row.names=1, threads=threads)
  taxonomy <- read.csv(taxonomy_file, row.names=1)
  metadata <- read.csv(metadata_file, row.names=1)
  tree <- read_tree(tree_file)

  metagenome <- phyloseq(
    otu_table(abundance, taxa_are_rows=FALSE),
    tax_table(as.matrix(taxonomy)),
    sample_data(metadata),
    tree
  )
  
  for (v in c('Site_id', 'Year', 'PCA_Grp', 'group')) {
    sample_data(metagenome)[[v]] <- factor(sample_data(metagenome)[[v]])
  }
  sample_data(metagenome)[['Eruption']] <- factor(sample_data(metagenome)[['Eruption']],
                                                  c('PreEruption', 'PostEruption'))
  sample_data(metagenome)[['Month']] <- factor(sample_data(metagenome)[['Month']],
                                                  c('March', 'August', 'November'))

  return(metagenome)
}

##===================== Metagenome filtering =====================##

filter_prevalence <- function(metagenome, min_prevalence=0) {
  otus <- taxa_names(metagenome)
  prevalence <- colSums(otu_table(metagenome) > 0)
  otus_to_keep <- otus[prevalence >= min_prevalence]
  
  metagenome <- prune_taxa(otus_to_keep, metagenome)
  metagenome <- prune_samples(sample_sums(metagenome)>0, metagenome)
  
  return(metagenome)
}

filter_samples <- function(metagenome, min_sample_size=0) {

  samples <- sample_names(metagenome)
  sample_sizes <- sample_sums(metagenome)
  samples_to_keep <- samples[sample_sizes >= min_sample_size]
  
  metagenome <- prune_samples(samples_to_keep, metagenome)
  metagenome <- prune_taxa(taxa_sums(metagenome)>0, metagenome)
  
  return(metagenome)
}

subset_metagenome <- function(metagenome, file, column='Genus') {
  otus <- taxa_names(metagenome)
  taxa <- read.table(file, header=TRUE)[, column]
  in_clade <- tax_table(metagenome)[, column] %in% taxa
  otu_clade <- otus[in_clade]
  metagenome <- prune_taxa(otu_clade, metagenome)
  
  return(metagenome)
}

filter_metagenome <- function(metagenome, min_sample_size=0, min_prevalence=0, 
                              taxa=NULL, clade_file=NULL) {

  metagenome <- subset_metagenome(metagenome, clade_file)
  metagenome <- filter_prevalence(metagenome, min_prevalence=min_prevalence)
  metagenome <- filter_samples(metagenome, min_sample_size=min_sample_size)
  
  if (!is.null(taxa)) {
    metagenome <- tax_glom(metagenome, taxa)
  }
  
  return(metagenome)
}

##======================== Misc functions ========================##

remove_outliers <- function(df) {
  for (column in colnames(df)) {
    outliers <- names(boxplot.stats(as.matrix(df)[, column])$out)
    df[outliers, column] <- NA
  }
  return(df)
}
