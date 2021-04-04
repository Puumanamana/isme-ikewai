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

group_colors <- c('01'='yellow', '02-04'='blue', '05-09'='green', '10-15'='purple', '16-20'='red')
eruption_colors <- c('salmon', 'darkcyan')
names(eruption_colors) <- c('PreEruption', 'PostEruption')

colors <- list(
  group=group_colors,
  Eruption=eruption_colors
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

