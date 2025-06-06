library <- function(...) suppressPackageStartupMessages(base::library(...))

library(dplyr)
library(tidyr)
library(Polychrome)
library(stringr)

rank <- snakemake@wildcards["rank"]
regex <- sprintf("%s__([^;]+)", substr(rank, 1, 1))

taxa <- lapply(snakemake@input, readLines) %>%
    unlist %>%
    str_match(regex) %>%
    `[`(,2) %>% unique %>% na.omit

palette <- createPalette(length(taxa), c("#00ff00", "#ff7800"))

colors <- data.frame(taxon = taxa, color = palette)
write.csv(colors, file = unlist(snakemake@output), row.names = F, quote = F)
