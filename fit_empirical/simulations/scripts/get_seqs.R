# Inputs:
# 1. Consensus sequence alignment including HXB2
# 2. Founder sequence
# 3. Output file name

# load library
library(ape)

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# read in consensus alignment and select sequences of interest
cons <- del.gaps(as.matrix(read.FASTA(args[1])[c('B.FR.83.HXB2_LAI_IIIB_BRU.K03455', 'CON_B(1295)')]))
#founder <- del.gaps(read.FASTA(args[2]))

# concatenate sequences
seqs <- cons
for(i in 2:(length(args)-1)){
  seqs <- c(seqs, del.gaps(read.FASTA(args[i])))
}

# write fasta
write.FASTA(seqs, args[length(args)])
