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
                            taxonomy_file='../data/taxonomy.csv',
                            metadata_file='../data/metadata.csv',
                            threads=4) {
  
  abundance <- fast_table_load(abundance_file, row.names=1, threads=threads)
  taxonomy <- read.csv(taxonomy_file, row.names=1)
  metadata <- read.csv(metadata_file, row.names=1)
  metadata$site_id <- factor(metadata$site_id)
  metadata$year <- factor(metadata$year)
  metadata$eruption <- factor(metadata$eruption, levels=c('pre-eruption', 'post-eruption'))
  
  metagenome <- phyloseq(
    otu_table(abundance, taxa_are_rows=FALSE),
    tax_table(as.matrix(taxonomy)),
    sample_data(metadata)
  )
  
  for (v in c('site_id', 'year')) {
    sample_data(metagenome)[[v]] <- factor(sample_data(metagenome)[[v]])
  }

  return(metagenome)
}

##===================== Metagenome filtering =====================##

filter_metagenome <- function(metagenome, min_prevalence=0, taxa=NULL, clade_file=NULL) {
  samples <- sample_names(metagenome)
  otus <- taxa_names(metagenome)
  
  prevalence <- colSums(otu_table(metagenome) > 0)
  otus_to_keep <- otus[prevalence >= min_prevalence]
  
  if (!is.null(clade_file)) {
    geni <- read.table(clade_file, header=TRUE)$Genus
    in_clade <- tax_table(metagenome)[,'Genus'] %in% geni
    otu_clade <- otus[in_clade]
    otus_to_keep <- intersect(otus_to_keep, otu_clade)
  }
  
  metagenome <- prune_taxa(otus_to_keep, metagenome)
  metagenome <- prune_samples(samples[sample_sums(metagenome)>0], metagenome)
  
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
