library(lhs)
library(tidyverse)
# set.seed(20250203)
set.seed(20250211)

param_ranges <- list(generation_time = c(1, 2),
                     mutation_rate = c(1e-5, 1e-4),
                     recombination_rate = c(1e-6, 1e-4),
                     conserved_cost = c(0, 0.99),
                     replicative_cost = c(0, 0.003),
                     immune_cost = c(0, 0.9),
                     immune_response_count = c(1, 1000),
                     days_full_potency = c(1, 1000),
                     active_to_latent = c(0.0001, 0.01), # latent to active c(0.001, 0.1)
                     latent_proliferate_die = c(0.001, 0.1), 
                     starting_pop_size = c(1, 1000), 
                     max_growth_rate = c(0.01, 1),
                     effective_pop_size = c(1000, 10000),
                     epitope_length = c(1, 20)) # change to 5-15?

# transform value in [0,1] range to value in other range
# x is in range [0,1]. want to transform back to true parameter value
# param_range is [range_min, range_max]
linearly_transform <- function(x, param_range){
  range_min <- param_range[1]
  range_max <- param_range[2]
  # (x-0)/(1-0)*(range_max - range_min) + range_min
  x*(range_max - range_min) + range_min
}

# sapply(param_ranges, function(x) sapply(seq(0, 1, 0.1), function(y) linearly_transform(y, x))) %>% View()

# subtype - don't have a good delta equivalent, hopefully won't be that sensitive in local analysis
# q_matrix  - don't have a good delta equivalent, and less sensitive than other parameters
# epitope_location - don't have a good delta equivalent, and less sensitive than other parameters
# number_of_epitopes - changes the immune cost values too; similar to epitope length for local sensitivity analysis
# reference sequence
# nubmer of conserved sites and location
# immune start day




# number of parameters
n_param <- length(param_ranges)
# number of trajectories
n_traj <- 3*n_param
# # number of replicates
# n_rep <- 1
# # number of simulations
# n_sim <- n_traj*(n_rep*(n_param+1))
# latin hypercube sampling 
# "optimum" seeks to maximize mean distance from each point to all other points
lhs <- optimumLHS(n_traj, n_param, maxSweeps = 2, eps = 0.1, verbose = TRUE)
colnames(lhs) <- names(param_ranges)
rownames(lhs) <- 1:n_traj

linearly_transform(lhs, c(-1,1)) 

lhs %>% as.data.frame() %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(x = value)) +
  facet_wrap(~name) +
  geom_histogram(bins = 30)

linearly_transform(lhs, c(-1,1)) %>% as.data.frame() %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(x = value)) +
  facet_wrap(~name) +
  geom_histogram(bins = 30)

2^linearly_transform(lhs, c(-1,1)) %>% as.data.frame() %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(x = value)) +
  facet_wrap(~name) +
  geom_histogram(bins = 30) +
  scale_x_continuous(trans = 'log2')

2^linearly_transform(lhs, c(-1,1))

base_values <- list(generation_time = 1,
                     mutation_rate = 4e-5,
                     recombination_rate = 2e-5,
                     conserved_cost = 0.5,
                     replicative_cost = 0.001,
                     immune_cost = 0.3,
                     immune_response_count = 100,
                     days_full_potency = 500,
                     active_to_latent = 0.001, # latent to active c(0.001, 0.1)
                     latent_proliferate_die = 0.001, 
                     starting_pop_size = 100, 
                     max_growth_rate = 0.05,
                     effective_pop_size = 4000,
                     epitope_length = 10) 




# delta <- sort(rep(c(0.1, 0.2, 0.3), nrow(lhs)/3))
delta <- 0.2 #rep(0.2, nrow(lhs))
trajectories <- lapply(1:nrow(lhs), function(row){
  lhs_row <- lhs[row,]
  traj <- rbind(c(lhs[row,], which_diff=NA, delta=NA), t(sapply(sample(1:length(lhs_row)), function(x){
    d <- sample(c(-delta, delta), 1) #sample(c(-delta[row], delta[row]), 1)
    old_val <- lhs_row[x]
    new_val <- old_val+d
    if(new_val < 0 | new_val > 1){
      new_val <- old_val-d
    }
    lhs_row[x] <<- new_val
    return(c(lhs_row, which_diff=colnames(lhs)[x], delta=d))
  })))
  rownames(traj) <- paste0(row, '_', 0:n_param)
  data.frame(traj) %>% 
    rownames_to_column() %>% 
    as_tibble() #%>% 
    # mutate(delta = delta[row])
}) %>% bind_rows() %>% 
  separate(rowname, into = c('traj', 'i'))



parameter_sets <- lapply(1:ncol(trajectories), function(x){
  transformed <- unlist(trajectories[,x])
  if(colnames(trajectories)[x] %in% names(param_ranges)){
    transformed <- linearly_transform(as.numeric(transformed), param_ranges[[x-2]])
    if(names(param_ranges)[x-2] %in% c('days_full_potency', 'starting_pop_size', 'effective_pop_size', 'epitope_length', 'immune_response_count')){
      transformed <- round(transformed)
    }
  }
  transformed <- tibble(transformed)
  names(transformed) <- colnames(trajectories)[x]
  return(transformed)
}) %>% bind_cols()
parameter_sets

write_csv(parameter_sets, '~/Documents/research/wavess_manuscript/morris_sensitivity/parameter_sets.csv')
