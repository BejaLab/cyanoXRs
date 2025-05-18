library <- function(...) suppressPackageStartupMessages(base::library(...))

library(dplyr)
library(tidyr)
library(treeio)
library(ggtree)
library(readxl)
library(ggstar)
library(ggplot2)
library(ggnewscale)
library(castor)
library(ape)
library(tools)
library(tibble)
library(ggtreeExtra)
library(tools)
library(adephylo)
library(tools)
library(readxl)
library(ggnewscale)
library(phangorn)
library(randomcoloR)

with(snakemake@input, {
    scaffolds_file <<- scaffolds
    tree_file <<- tree
    synonyms_file <<- synonyms
    ref_taxonomy_file <<- ref_taxonomy
})
output_file <- unlist(snakemake@output)
with(snakemake@params, {
    root_by <<- root_by
})

extract_tax <- function(.data, col, taxon, prefix) {
    extract(.data, col, into = taxon, regex = sprintf("%s__(.+?);", prefix), remove = F)
}

synonyms <- read.table(synonyms_file, col.names = c("from", "to")) %>%
    with(setNames(to, from))
scaffolds <- data.frame(line = readLines(scaffolds_file)) %>%
    separate(line, into = c("scaffold", "scaf_description"), sep = " ", extra = "merge") %>%
    mutate(scaffold = sub(":", "_", scaffold))
refs <- read.table(ref_taxonomy_file, col.names = c("label", "ref_taxonomy"))

to_treedata <- function(tree) {
    class(tree) <- c("tbl_tree", "tbl_df", "tbl", "data.frame")
    as.treedata(tree)
}

tree <- read.tree(tree_file) %>%
    as_tibble %>%
    mutate(is.tip = ! node %in% parent) %>%
    mutate(scaffold = sub("_\\d+$", "", label)) %>%
    left_join(scaffolds, by = "scaffold") %>%
    left_join(refs, by = "label") %>%
    mutate(description = paste(scaf_description, ref_taxonomy)) %>%
    extract_tax("description", "domain", "d") %>%
    extract_tax("description", "phylum", "p") %>%
    mutate(phylum = recode(phylum, !!!synonyms)) %>%
    mutate(label.show = ifelse(label %in% refs$label, label, NA)) %>%
    mutate(taxon = ifelse(is.na(phylum), domain, phylum)) %>%
    group_by(taxon) %>%
    mutate(taxon = ifelse(n() > 1, taxon, NA)) %>%
    separate(label, into = c("SH_aLRT", "UFboot"), sep = "/", fill = "left", remove = F) %>%
    mutate(UFboot = ifelse(is.tip | UFboot == "", NA, UFboot)) %>%
    mutate(UFboot = as.numeric(UFboot))

n_taxa <- filter(tree, !is.na(taxon)) %>%
    distinct(taxon) %>% nrow
palette <- distinctColorPalette(n_taxa)

p <- ggtree(to_treedata(tree), aes(color = taxon), layout = "circular", linewidth = 0.4) +
    geom_tiplab(aes(label = label.show), size = 4, offset = 0.1, align = F) +
    geom_point2(aes(subset = !is.na(UFboot) & UFboot >= 95), size = 1, color = "darkgray") +
    scale_color_manual(values = palette) +
    geom_treescale(width = 0.5)

ggsave(output_file, p, width = 12, height = 12)
