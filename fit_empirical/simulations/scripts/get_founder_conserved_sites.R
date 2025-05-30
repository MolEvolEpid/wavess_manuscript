# Inputs:
# 1. HXB2 conserved sites
# 2. Mapped HXB2-founder positions (hxb2_founder_positions.tsv)
# 3. Output file name for founder conserved sites
# 4. Output file name for founder conserved sites plot

# load library
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ape))
library(wavess)

flt_aln <- read.FASTA(snakemake@input[[1]])
ref_founder_aln <- read.FASTA(snakemake@input[[2]])

identify_conserved_sites(flt_aln, 'sim_founder', ref = 'B.FR.83.HXB2_LAI_IIIB_BRU.K03455', founder_aln = ref_founder_aln) |>
  filter(conserved == "Yes") |>
  select(founder_pos, founder_base) |>
  deframe() |>
  toupper() |>
  enframe(name = "position", value = "nucleotide") |>
  write_csv(snakemake@output[[1]])

