library <- function(...) suppressPackageStartupMessages(base::library(...))

library(dplyr)
library(tidyr)
library(treeio)
library(ggtree)
library(ggplot2)
library(readxl)
library(castor)
library(ape)
library(tibble)

with(snakemake@input, {
    tree_file <<- tree
    metadata_file <<- metadata
    refs_file <<- refs
})
with(snakemake@output, {
    image_file <<- image
    jtree_file <<- jtree
})

with(snakemake@params, {
    ingroup <<- ingroup
})

to_treedata <- function(tree) {
    class(tree) <- c("tbl_tree", "tbl_df", "tbl", "data.frame")
    as.treedata(tree)
}

refs <- names(read.fasta(refs_file)) %>%
    sub(" .+", "", .)
outgroups <- refs[! refs %in% ingroup]

metadata <- read_excel(metadata_file) %>%
    rename(taxonomy = `GTDBtk taxonomy`) %>%
    extract(taxonomy, into = c("Family", "Genus"), regex = "(o__[^;]+;f__[^;]*);(g__[^;]*)") %>%
    mutate(label_cleaned = sub(":", "_", label))

tree <- read.tree(tree_file) %>%
    root(outgroups[outgroups %in% .$tip.label], edgelabel = T, resolve.root = T) %>%
    drop.tip(outgroups) %>%
    as_tibble %>%
    mutate(is.tip = ! node %in% parent & node != parent) %>%
    mutate(Support = ifelse(is.tip, NA, as.numeric(label))) %>%
    rename(label_cleaned = label) %>%
    left_join(metadata, by = "label_cleaned") %>%
    mutate(Label.show = case_when(
        type == "isolate" ~ strain,
        type == "MAG" & !is.na(Genus) & Genus != "g__" ~ sprintf("%s (%s)", assembly, Genus),
        type == "MAG" & !is.na(Family) ~ assembly,
        F ~ NA)
    ) %>%
    to_treedata

shapes <- c(
    freshwater = 24, # triangle up
    marine = 25, # trinagle down
    terrestrial = 22, # square
    saline = 23 # diamond
)
colors <- c(
    cold = "#55ddff",
    hot = "#ff5555",
    moderate = "#2aff80"
)

gap <- 0.03
p <- ggtree(tree, aes(color = Family), size = 0.1, layout = "rectangular", open.angle = 15) +
    geom_tiplab(aes(label = Label.show), size = 2, align = F) +
    geom_text2(aes(subset = Support > 60, x = branch, label = Support), color = "black", size = 2, vjust = -0.5) +
    geom_tippoint(aes(subset = !is.na(habitat1), shape = habitat1, fill = habitat2, x = x + gap * 1), color = "transparent") +
    scale_shape_manual(values = shapes, na.value = 21) +
    scale_fill_manual(values = colors) +
    geom_treescale(width = 0.1)

write.jtree(tree, jtree_file)
ggsave(image_file, p, width = 9, height = 7.5)
