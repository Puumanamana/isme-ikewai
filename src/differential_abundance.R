
#============ Set up environment ============#

library(phyloseq)
library(DESeq2)
library(vegan)
library(BiocParallel)
library(ggplot2)

source('loader.R')

#==================================================#
#--------------------- CLI ------------------------#
#==================================================#
library(argparse)

parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option
parser$add_argument("--rds", help='metagenome in RDS format', default='data/raw_metagenome.RDS')
parser$add_argument("--h5", help='metagenome in h5 format', default='data/h5')
parser$add_argument("-o", "--outdir", help='output folder', default='results')
parser$add_argument("--clade-dir", default='data/functional_clades', help='path to directory with clade files')
parser$add_argument("-c", "--clade", default='all')
parser$add_argument("-f", "--factors", default=c('PCA_Chem_Group', 'eruption'), nargs='+')
parser$add_argument("-t", "--threads", type='integer', default=2)
parser$add_argument("--pval", type='double', default=0.05)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

clade_file <- NULL
if (args$clade %in% c('nirS', 'dsrB')) {
  clade_file <- sprintf('%s/%s_genuses.tsv', args$clade_dir, args$clade)
}

factors_str <- paste(args$factors, collapse='-')

register(MulticoreParam(args$threads))

#=========== Loading and Filtering ===========#
metagenome <- load_h5(args$h5)
metagenome <- subset_samples(metagenome, PCA_Chem_Group %in% c('C', 'D'))
metagenome <- filter_data(metagenome, min_prevalence=1, clade_file=clade_file)

#==================================================#
#---------------- Simper (vegan) ------------------#
#==================================================#

abundance <- as.data.frame(otu_table(metagenome))
factor_values <- sample_data(metagenome)[[args$factors[1]]]

strata <- NULL
if (length(args$factors) > 1) {
  strata <- sample_data(metagenome)[[args$factors[2]]]
}
simper_obj <- simper(abundance, factor_values, permutations=99, parallel=args$threads,
                     strata=strata)

simper_out <- sprintf("%s/differential-otus/simper/%s/%s", 
                      args$outdir, args$clade, factors_str)
dir.create(simper_out, recursive=TRUE, showWarnings=FALSE)

for (comparison in names(simper_obj)) {
  simper_table <- as.data.frame(simper_obj[[comparison]])
  write.csv(simper_table, quote=FALSE,
            file=sprintf('table-%s/%s.csv', simper_out, comparison))
  simper_table <- simper_table[simper_table$p < args$pval, ]
  write.csv(simper_table, quote=FALSE,
            file=sprintf('list-%s/%s.csv', simper_out, comparison))
}

#==================================================#
#---------------- DESeq 2 package -----------------#
#==================================================#

## Script from tutorial at http://joey711.github.io/phyloseq-extensions/DESeq2.html

#============= Fit DESeq2 model ===============#
design <- as.formula(sprintf(" ~ %s", paste(args$factors, collapse='+')))

diagdds <- phyloseq_to_deseq2(metagenome, design)

# calculate geometric means prior to estimate size factors
# DESeq2 cannot normalize when every OTU contains at least one zero
geom_mean <- function(x) { exp(mean(log(x[x>0]), na.rm=TRUE)) }

geoMeans <- apply(counts(diagdds), 1, geom_mean)
diagdds <- estimateSizeFactors(diagdds, geoMeans=geoMeans)
diagdds <- DESeq(diagdds, fitType="local", test="Wald", parallel=TRUE, quiet=TRUE)

#========== Filter and sort result ============#

sig_otu_table <- results(diagdds, alpha=args$pval, pAdjustMethod='BH',parallel=TRUE)
sig_otu_table <- sig_otu_table[order(sig_otu_table$padj, na.last=NA), ]

sig_otu_table <- cbind(
  as(sig_otu_table, "data.frame"), 
  as(tax_table(metagenome)[rownames(sig_otu_table), ], "matrix")
)

deseq2_out <- sprintf("%s/differential-otus/deseq2", args$outdir)
dir.create(sprintf("%s/table", deseq2_out), recursive=TRUE, showWarnings=FALSE)
dir.create(sprintf("%s/list", deseq2_out), recursive=TRUE, showWarnings=FALSE)

write.csv(sig_otu_table, quote=FALSE,
          file=sprintf('%s/table/%s_by_%s.csv', deseq2_out, args$clade, factors_str))

sig_otus <- rownames(sig_otu_table[order(-sig_otu_table$padj), ])
write.table(sig_otus, quote=F, col.names="OTU", row.names=F,
            file=sprintf('%s/list/%s_by_%s.csv', deseq2_out, args$clade, factors_str))

#==================================================#
#------------ Plot OTU distributions --------------#
#==================================================#

abundance_sub <- as.data.frame(otu_table(metagenome)[, differential_otus])

for (fact in args$factors) {
  abundance_sub[[fact]] <- sample_data(metagenome)[[fact]]
}
data <- reshape2::melt(abundance_sub, id_vars=args$factors)

fill <- data[[args$factors[1]]]
if (length(args$factors) == 1) {
  x <- fill
} else {
  x <- data[[args$factors[2]]]
}

ggplot(data, aes(x=x, y=value, fill=fill)) + 
  geom_boxplot() +
  scale_y_continuous(trans=scales::pseudo_log_trans(base=10)) +
  facet_wrap('~ variable', scale="free") +
  labs(x='', y='abundance', fill=args$factors[1])

ggsave(sprintf('%s/distr_%s_%s.pdf', args$outdir, factors_str, args$clade), 
       width=40, height=30, units="cm")

