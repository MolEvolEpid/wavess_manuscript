# Load libraries
suppressPackageStartupMessages(library(tidyverse))
library(wavess)
#library(ggtree)
suppressPackageStartupMessages(library(ape))
library(phytools)

# Read in data
founder_len <- ncol(as.matrix(read.FASTA(snakemake@input[[1]])))
tr <- read.tree(snakemake@input[[2]])
#tr <- unroot(drop.tip(tr, 'CON_B_1295_'))
first_gen <- min(as.numeric(gsub('gen|_.*', '', tr$tip.label)))
earliest_samps <- tr$tip.label[grep(paste0('gen', first_gen, '_'), tr$tip.label)]

if(length(earliest_samps) == 1){
tr <- tr %>%
    reroot(which(tr$tip.label == earliest_samps))
}else{
tr <- tr %>%
    reroot(getMRCA(tr, earliest_samps))
}

#gens <- gsub("gen|_.*", "", labels(aln))
#names(gens) <- labels(aln)
#calc_div_metrics(aln, "founder0", gens) |>
#  mutate(file = snakemake@output[[1]]) |>
#  write_csv(snakemake@output[[1]])

gens <- gsub("gen|_.*", "", tr$tip.label)
names(gens) <- tr$tip.label
gens <- factor(gens, levels = sort(unique(as.numeric(gens))))

calc_tr_stats(tr, gens, 1/founder_len) |>
  mutate(file = snakemake@output[[1]]) |>
  write_csv(snakemake@output[[1]])

#cols <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000', 'green', 'red', 'blue', 'orange', 'purple', 'yellow', 'cyan', 'grey', 'pink', 'cadetblue', 'forestgreen', 'yellow', 'darkseagreen', 'darkslategray')

#tr_plot <- ggtree(tr) +
#  geom_tippoint(aes(col = factor(c(gens, rep(NA, Nnode(tr))), levels = sort(unique(as.numeric(gens)))))) +
#  scale_color_manual(values = cols) +
#  geom_treescale() +
#  labs(col = "Generation")

#ggsave(snakemake@output[[2]], tr_plot)

