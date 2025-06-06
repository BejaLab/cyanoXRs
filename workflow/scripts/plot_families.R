library <- function(...) suppressPackageStartupMessages(base::library(...))

library(dplyr)
library(tidyr)
library(ComplexUpset)
library(scales)
library(ggplot2)

with(snakemake@input, {
    blast_file <<- blast
    scaffolds_file <<- scaffolds
    hmmsearch_file <<- hmmsearch
    taxa_file <<- taxa
})
with(snakemake@params, {
    min_id <<- min_id
    max_evalue <<- evalue
})
with(snakemake@output, {
    counts_file <<- counts
    plot_file <<- plot
})

taxa <- read.csv(taxa_file)
get_taxon <- function(.data, rank) {
    regex <- sprintf("%s__([^;]+);", substr(rank, 1, 1))
    extract(.data, genome, into = rank, regex = regex, remove = F)
}
scaffolds <- readLines(scaffolds_file) %>%
    data.frame %>%
    separate(1, into = c("scaffold", "genome"), sep = " ", extra = "merge") %>%
    mutate(genome = sub("megahit:[^ ]+", "", genome)) %>%
    get_taxon("order") %>% get_taxon("family")

hmmsearch <- readLines(hmmsearch_file) %>%
    data.frame %>%
    extract(1, into = "qseqid", regex = "^(\\S+)") %>%
    filter(qseqid != "#")

blast_cols <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle")
matches <- read.table(blast_file, col.names = blast_cols, sep = "\t") %>%
    filter(pident >= min_id, evalue <= max_evalue) %>%
    bind_rows(hmmsearch) %>%
    distinct(qseqid, .keep_all = T) %>%
    mutate(scaffold = sub("_[0-9]+$", "", qseqid)) %>%
    extract(stitle, into = "clade", regex = "/([^/]+)/") %>%
    replace_na(list(clade = "unclassified"))

genomes <- left_join(scaffolds, matches, by = "scaffold") %>%
    group_by(genome) %>%
    filter(!is.na(clade) | all(is.na(clade))) %>%
    distinct(genome, clade, order)
write.csv(genomes, file = counts_file, na = "absent")

presence <- mutate(genomes, present = 1) %>%
    group_by(order) %>%
    mutate(order_label = sprintf("%s (%s/%s)", order, n_distinct(genome[!is.na(clade)]), n_distinct(genome))) %>%
    left_join(taxa, by = c(order = "taxon")) %>%
    spread(clade, present, fill = 0)

colors <- distinct(presence, order_label, color) %>%
    filter(!is.na(color)) %>%
    with(setNames(color, order_label))
clades <- ungroup(genomes) %>%
    filter(!is.na(clade)) %>%
    distinct(clade) %>% pull
orders_p <- ggplot(map = aes(fill = order_label)) +
    geom_bar(stat = 'count', position = 'fill') +
    scale_y_continuous(labels = percent_format()) +
    scale_fill_manual(values = colors) +
    ylab('Orders')
p <- upset(presence, clades, annotations = list(Orders = orders_p))

ggsave(plot_file, p)
