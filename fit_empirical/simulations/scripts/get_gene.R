# Inputs:
# 1. Alignment including hxb2, consensus, and founder
# 2. HXB2 start position of gene
# 3. HXB2 end position of gene
# 4. Output file name

# load library
library(ape)

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# sequence names (to make more generalizable, could include as input arguments)
hxb2 <- 'B.FR.83.HXB2_LAI_IIIB_BRU.K03455'
consensus <- 'CON_B(1295)'
founder <- 'B.US.2011.DEMB11US006.KC473833'

# read in alignments
aln <- as.matrix(read.FASTA(args[1]))

n_not_gap <- cumsum(as.character(aln[hxb2,]) != '-')
start <- which(n_not_gap == args[2])
end <- which(n_not_gap == args[3])

aln_sub <- aln[,start:end]

if(nrow(aln_sub) == 4){
  founder_orig <- founder
  founder <- labels(aln_sub)[4]
  cum_not_dash <- apply(as.character(aln_sub) != '-', 1, cumsum)[,4]
  founder_first_pos <- which(cum_not_dash == 1)
  founder_last_pos <- min(which(cum_not_dash == max(cum_not_dash)))
  sim_founder_seq <- as.matrix(as.DNAbin(c(as.character(aln_sub[founder_orig,][1:(founder_first_pos-1)]),
		   as.character(aln_sub[founder,][founder_first_pos:founder_last_pos]),
		   as.character(aln_sub[founder_orig,][(founder_last_pos+1):ncol(aln_sub)]))))

  if(founder_first_pos == 1){
    sim_founder_seq <- as.matrix(as.DNAbin(c(
                   as.character(aln_sub[founder,][founder_first_pos:founder_last_pos]),
                   as.character(aln_sub[founder_orig,][(founder_last_pos+1):ncol(aln_sub)]))))
  }

  founder_real <- founder
  founder <- 'sim_founder'
  rownames(sim_founder_seq) <- founder
  aln_sub <- rbind(aln_sub, sim_founder_seq)
}

aln_no_gaps_in_founder <- aln_sub[,as.character(aln_sub[founder,]) != '-']

write.FASTA(aln_sub, paste0(args[4], '_all_aligned.fasta'))
write.FASTA(del.colgapsonly(aln_sub[c(hxb2,founder),]), paste0(args[4], '_hxb2_founder_aligned.fasta'))
rep(list(as.character(aln_no_gaps_in_founder[founder,])[1,]), 10) |>
  as.DNAbin() |>
  as.matrix() |>
write.FASTA(paste0(args[4], '_founder.fasta'))
write.FASTA(aln_no_gaps_in_founder[consensus,], paste0(args[4], '_cons_no_founder_gaps.fasta'))

