library(wavess)
library(ape)
library(readr)

epi_probs <- get_epitope_frequencies(env_features$Position)

ref_founder_map <- map_ref_founder(read.FASTA(snakemake@input[[1]]),
  ref = "B.FR.83.HXB2_LAI_IIIB_BRU.K03455",
  founder = "sim_founder")

sample_epitopes(epi_probs, max_fit_cost = as.numeric(snakemake@wildcards$max_epi), ref_founder_map = ref_founder_map) |>
  write_csv(snakemake@output[[1]])
