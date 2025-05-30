# Input:
# 1. Fasta file to separate into individual files
# 2. Output file name

# load library
library(ape)

# read in consensus alignment and select sequences of interest
seqs <- read.FASTA(snakemake@input[[1]])

# save each sequence to a separate file
for(seq_name in names(seqs)){
  name_info <- strsplit(seq_name, split = '\\+')[[1]]
  name_new <- paste0(name_info[1], "+1+B")
  start_pos <- as.numeric(name_info[2])
  if(name_new == "105692+1+B" & start_pos > 6000) name_new <- "105692+197+B"
  if(name_new == "105693+1+B" & start_pos > 6000) name_new <- "105693+202+B"
  if(name_new == "105695+1+B" & start_pos > 6000) name_new <- "105695+181+B"
  if(name_new == "105696+1+B" & start_pos > 6000) name_new <- "105696+207+B"
  if(name_new == "105698+1+B" & start_pos > 6000) name_new <- "105698+280+B"
  dir.create(paste0("output/", name_new, "/inputs/"), showWarnings = FALSE, recursive = TRUE)
  write.FASTA(seqs[seq_name], paste0("output/", name_new, "/inputs/", name_new, "_original_founder.fasta"))
}

