library <- function(...) suppressPackageStartupMessages(base::library(...))

library(treeio)
library(ggtree)
library(dplyr)
library(tidyr)
library(stringr)
library(ape)
library(phangorn)

with(snakemake@input, {
    genomes_file <<- genomes
    taxonomy_file <<- taxonomy
    tree_file <<- tree
    colors_file <<- colors
})
output_file <- unlist(snakemake@output)

rank_prefix <- "o__"
rank_regex <- paste0(rank_prefix, "[^;]+")

tree <- read.tree(tree_file)
tree$node.label <- gsub("'", "", str_extract(tree$node.label, rank_regex))

colors <- read.csv(colors_file)

taxa <- read.csv(genomes_file) %>%
    separate(genome, into = c("assembly", "lineage"), sep = " ", extra = "merge") %>%
    group_by(order) %>%
    summarize(n_assemblies = n_distinct(assembly), n_absent = sum(clade == "absent"), assembly = paste(unique(assembly), collapse = ";")) %>%
    left_join(colors, by = c(order = "taxon")) %>%
    mutate(taxon = paste0(rank_prefix, order))

taxonomy <- read.table(taxonomy_file, col.names = c("assembly", "lineage"), sep = "\t") %>%
    filter(assembly %in% tree$tip.label) %>%
    mutate(taxon = str_extract(lineage, rank_regex)) %>%
    group_by(taxon) %>%
    summarize(n = n(), assembly = paste(assembly, collapse = ";"))

singletons <- filter(taxonomy, n == 1) %>%
    select(taxon, assembly)

taxa_labels <- with(tree, node.label[!is.na(node.label)])

for (tax in taxa_labels) {
    node <- with(tree, length(tip.label) + which(node.label == tax))
    tips <- unlist(Descendants(tree, node))
    dists <- dist.nodes(tree)[node,tips]
    dists <- dists[order(dists)]
    tips <- names(dists[-c(1, length(tips))])
    tree <- drop.tip(tree, as.numeric(tips), collapse.singles = F)
}

singleton_list <- with(singletons, setNames(assembly, taxon))
for (tax in names(singleton_list)) {
    assembly <- singleton_list[tax]
    node <- with(tree, which(tip.label == assembly))
    dummy <- read.tree(text = sprintf("(%s_1:0,%s_2:0)%s:0;", assembly, assembly, tax))
    tree <- bind.tree(tree, dummy, where = node)
}

p <- ggtree(tree)
p$data <- left_join(p$data, singletons, by = c(label = "assembly")) %>%
    mutate(taxon = ifelse(is.na(taxon), label, taxon)) %>%
    left_join(taxa, by = "taxon") %>%
    mutate(label = ifelse(is.na(taxon), "", sprintf("%s (%d/%d)", taxon, n_assemblies - n_absent, n_assemblies)))

chosen_taxa <- filter(p$data, grepl(rank_prefix, label), !isTip) %>%
    split(f = 1:nrow(.))

for (tax in chosen_taxa) {
    p <- collapse(p, node = tax$node, "min", fill = tax$color)
}

colors <- filter(p$data, grepl(rank_prefix, label)) %>%
    with(setNames(color, label))

p <- p +
    geom_nodelab(hjust = 0) +
    geom_tiplab() +
    geom_point2(aes(subset = grepl(rank_prefix, label), color = label), shape = 15, size = 5) +
    scale_color_manual(values = colors) +
    geom_treescale(width = 0.1) +
    theme(legend.position = "none")

ggsave(output_file, p, width = 4, height = 6)
