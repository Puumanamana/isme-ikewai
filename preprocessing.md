# Preprocessing for each plot

### Global preprocessing steps

Discard samples with no metadata information
Discard OTUs with a null abundance (of course)

### Figure 5 (CCA pre- and post-eruption) and Figure S3 (CCA per PCA Chem Group)

1. Discard OTUs that appear in less than 2 samples
2. Subsample at 5,000 reads per sample, discard samples with less than this amount
3. Convert to relative abundance

### Figure 6 (Alpha diversity)

1. Discard OTUs that appear in less than 2 samples

### Figure 7 (Taxa barplots) and S4 (Taxa barplots for dsrA/nirS by (site_id, eruption))

1. Subset PCA Chem groups C and D only
2. Aggregate OTU per genus
3. Sum counts per factor level (e.g. PCA_Chem_Group or eruption) and normalize to frequencies

### Figure 8 (Boxplot and DESeq2 scores)

1. Subset PCA Chem groups C and D only
2. Aggregate OTU per genus
3. Discard OTUs that appear in less than 2 samples

