# Load libraries
suppressPackageStartupMessages(library(tidyverse))
library(ape)

# Read in data
seqs <- read.dna(snakemake@input[[1]], format = 'fasta', as.matrix = TRUE)
aln <- as.matrix(read.FASTA(snakemake@input[[2]]))
cons_seq <- as.matrix(read.FASTA(snakemake@input[[3]]))

# Find part of sequence corresponding to real data
cum_not_dash <- apply(as.character(aln) != '-', 1, cumsum)[,4]
founder_first_pos <- which(cum_not_dash == 1)
placeholder_gaps <- sum(as.character(aln[3,1:founder_first_pos]) == '-')
founder_first_pos <- founder_first_pos - placeholder_gaps
founder_len <- sum(as.character(aln[4,]) != '-')
founder_last_pos <- founder_first_pos + founder_len - 1

# Subset to real person part of sequence
cons_seq <- cons_seq[,founder_first_pos:founder_last_pos]
seqs <- seqs[,founder_first_pos:founder_last_pos]

# Save fasta with no founder or consensus (for tree building)
write.FASTA(seqs[!grepl('founder', rownames(seqs)),], snakemake@output[[1]])

# Save fasta with founder but no consensus (for summary stats)
#write.FASTA(seqs, snakemake@output[[2]])

# Save fasta with consensus but no founder (for alternate tree building)
#write.FASTA(rbind(cons_seq, seqs[!grepl('founder', rownames(seqs)),]), snakemake@output[[3]])
