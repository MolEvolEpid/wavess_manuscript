library(tidyverse)
library(ape)

seqs_to_keep <- read_csv(snakemake@input[[1]]) %>%
      mutate(timepoints = gsub('^-', 'neg', timepoints)) %>%
      separate_rows(c(timepoints, seq_in_timepoint), sep = '-') %>%
      separate_rows(seq_in_timepoint, sep = '\\+') %>%
        mutate(timepoints = as.numeric(gsub('neg', '-', timepoints)),
               patient_id = gsub('\\+.*', '', aln_ids)) %>% 
      select(patient_id, aln_id = aln_ids, timepoint = timepoints, accession_number = seq_in_timepoint) %>%
  mutate(to_rm = case_when(patient_id == 10586 & timepoint > 4200 ~ 'rm',
                           patient_id == 13889 & timepoint > 4200 ~ 'rm',
                           TRUE ~ 'keep')) %>% 
  filter(to_rm == 'keep') %>% 
  pull(accession_number)

aln <- read.FASTA(snakemake@input[[2]])

accessions <- sapply(strsplit(labels(aln), '\\+'), function(y) y[4])
write.FASTA(del.gaps(aln[accessions %in% seqs_to_keep | is.na(accessions)]), snakemake@output[[1]]) # "Consensus_Subtype_B" becomes NA

