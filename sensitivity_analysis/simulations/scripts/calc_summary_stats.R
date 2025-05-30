# Load libraries
suppressPackageStartupMessages(library(tidyverse))
library(ape)
library(wavess)
library(ggtree)
library(phytools)

# Read in data
#aln <- as.matrix(read.FASTA(snakemake@input[[1]]))
#aln <- aln[grepl('founder0|gen', labels(aln)),]
tr <- read.tree(snakemake@input[[1]])
keepers <- tr$tip.label[tr$tip.label == 'founder0' | !grepl('founder', tr$tip.label)]
tr <- keep.tip(tr, keepers) |> reroot(which(tr$tip.label == "founder0")) |> drop.tip("founder0")

#gens <- gsub("gen|_.*", "", labels(aln))
#names(gens) <- labels(aln)
#calc_div_metrics(aln, "founder0", gens) |>
#  mutate(file = snakemake@output[[1]]) |>
#  write_csv(snakemake@output[[1]])

gens <- gsub("gen|_.*", "", tr$tip.label)
names(gens) <- tr$tip.label
calc_tr_stats(tr, gens, 1/1503) |>
  mutate(file = snakemake@output[[1]]) |>
  write_csv(snakemake@output[[1]])

#cols <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000', 'green', 'red', 'blue', 'orange', 'purple', 'yellow', 'cyan', 'grey', 'pink', 'cadetblue', 'forestgreen', 'yellow', 'darkseagreen', 'darkslategray')

#tr_plot <- ggtree(tr) +
#  geom_tippoint(aes(col = factor(c(gens, rep(NA, Nnode(tr))), levels = c('founder0', sort(unique(as.numeric(gens))))))) +
#  scale_color_manual(values = cols) +
#  geom_treescale() +
#  labs(col = "Generation")
#
#ggsave(snakemake@output[[3]], tr_plot)

