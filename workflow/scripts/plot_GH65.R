library <- function(...) suppressPackageStartupMessages(base::library(...))

library(treeio)
library(ggtree)
library(dplyr)
library(tidyr)
library(ape)

tree_file <- unlist(snakemake@input)
with(snakemake@params, {
    outgroup <<- outgroup
})

tree <- read.tree(tree_file)
tree$tip.label <- sub("_(\\d+)-(\\d+)_(.)_$", ":\\1-\\2(\\3)", tree$tip.label)
tree$tip.label <- sub("\\(_\\)$", "\\(+\\)", tree$tip.label)

tree <- ape::root(tree, outgroup, edgelabel = T, resolve.root = T)

p <- ggtree(tree) + geom_tiplab() + geom_text2(aes(x = branch, subset = as.numeric(label) >= 90, label = label)) + geom_treescale(width = 0.1)

ggsave(output_file, p, width = 3, height = 4)
