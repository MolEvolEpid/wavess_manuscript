import pandas as pd
import yaml

# args:
# config file
# nt_sub_probs
# output file path

df = pd.read_csv(snakemake.input[0])
row = df[(df['traj']==int(snakemake.wildcards.traj)) & (df['i']==int(snakemake.wildcards.i))]

with open(snakemake.input[1]) as istream:
    ymldoc = yaml.safe_load(istream)
    ymldoc['input_files']['inf_pop_size'] = snakemake.input[2]
    #ymldoc['input_files']['samp_scheme'] = snakemake.params['samp_scheme']
    ymldoc['input_files']['founder_seqs'] = snakemake.input[3]
    #ymldoc['input_files']['q'] = snakemake.params['q']
    #ymldoc['input_files']['conserved_sites'] = snakemake.params['conserved_sites']
    #ymldoc['input_files']['ref_seq'] = snakemake.params['ref_seq']
    ymldoc['input_files']['epitope_locations'] = snakemake.input[4]
    ymldoc['parameters']['generation_time'] = row['generation_time'].iloc[0].item()
    ymldoc['parameters']['mut_rate'] = row['mutation_rate'].iloc[0].item()
    ymldoc['parameters']['recomb_rate'] = row['recombination_rate'].iloc[0].item()
    ymldoc['parameters']['act_to_lat'] = row['active_to_latent'].iloc[0].item()
    ymldoc['parameters']['lat_to_act'] = row['active_to_latent'].iloc[0].item()*10
    ymldoc['parameters']['lat_prolif'] = row['latent_proliferate_die'].iloc[0].item()
    ymldoc['parameters']['lat_die'] = row['latent_proliferate_die'].iloc[0].item()
    ymldoc['parameters']['conserved_cost'] = row['conserved_cost'].iloc[0].item()
    ymldoc['parameters']['replicative_cost'] = row['replicative_cost'].iloc[0].item()
    #ymldoc['parameters']['immune_start_day'] = row['immune_start_day'].iloc[0].item()
    ymldoc['parameters']['n_for_imm'] = row['immune_response_count'].iloc[0].item()
    ymldoc['parameters']['time_to_full_potency'] = row['days_full_potency'].iloc[0].item()

with open(snakemake.output[0], 'w') as ostream:
    yaml.dump(ymldoc, ostream, default_flow_style=False, sort_keys=False)

