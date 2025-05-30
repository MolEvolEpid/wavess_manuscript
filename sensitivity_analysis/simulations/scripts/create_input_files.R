library(tidyverse)
library(wavess)
library(ape)

param_sets <- read_csv(snakemake@input[[1]])
founder <- as.matrix(read.FASTA(snakemake@input[[2]]))
#epitopes <- read_csv(snakemake@input[[3]])

epi_probs <- get_epitope_frequencies(env_features$Position)
gp120 <- slice_aln(ape::as.matrix.DNAbin(hxb2_cons_founder), 6225, 7787)
ref_founder_map <- map_ref_founder(gp120,
  ref = "B.FR.83.HXB2_LAI_IIIB_BRU.K03455",
  founder = "B.US.2011.DEMB11US006.KC473833"
)
epitopes <- sample_epitopes(epi_probs, aa_epitope_length = 15, ref_founder_map = ref_founder_map)
write_csv(epitopes, 'input/epitopes.csv')

define_sampling_scheme() |>
  write_csv('input/samp_scheme.csv')

for(x in 1:nrow(param_sets)){

  print(x)

  traj <- param_sets$traj[x]
  i <- param_sets$i[x]
  file_pref <- paste0('traj', traj, '_set', i)


  starting_pop_size <- param_sets$starting_pop_size[x]
  max_growth_rate <- param_sets$max_growth_rate[x]
  effective_pop_size <- param_sets$effective_pop_size[x]

  define_growth_curve(n0 = starting_pop_size, 
      max_growth_rate = max_growth_rate, 
      carry_cap = effective_pop_size, 
      n_gens = 10000) |>
    write_csv(paste0('input/growth_curves/', file_pref, '.csv'))


    new_founder <- founder
  for(n in 2:starting_pop_size){
    new_founder <- rbind(new_founder, founder)
  }
  write.FASTA(new_founder, paste0('input/founders/', file_pref, '.csv'))

  
  immune_cost <- param_sets$immune_cost[x]
  epitope_length <- param_sets$epitope_length[x]

  epitopes |>
    mutate(epi_end_nt = epi_start_nt+3*epitope_length,
	   max_fitness_cost = seq(0, immune_cost, length.out = 10+1)[2:(10+1)]) |>
    write_csv(paste0('input/epitopes/', file_pref, '.csv'))
}
