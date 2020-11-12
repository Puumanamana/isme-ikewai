library(ggplot2)
library(ggsci)
library(plyr)
library(dplyr)
library(tidyr)

io <- c(
  figures="../figures",
  nirS='../data/functional-clades/nirS_genuses.tsv',
  dsrB='../data/functional-clades/dsrB_genuses.tsv'
)

PCA_colors <- c(A='red', B='blue', C='green', D='purple')
eruption_colors <- c('salmon', 'darkcyan')
names(eruption_colors) <- c('pre-eruption', 'post-eruption')

colors <- list(
  PCA_Chem_Group=PCA_colors,
  eruption=eruption_colors
)

gginit <- function() {
  return(
    theme_linedraw() + 
    theme(
        title=element_text(size=16),
        strip.text.x=element_text(size=14),
        panel.grid=element_line("white"),
        axis.line=element_line("gray25"),
        axis.text=element_text(size=12, color="gray25"),
        axis.title=element_text(size=12, color="gray25"),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12)
    )
  )
}

