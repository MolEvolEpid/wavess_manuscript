from sys import argv
import yaml

# args:
# 1. config file
# 2. generation_info path
# 3. founder sequences
# 4. epitopes path
# 5. conserved sites path
# 6. reference sequence path
# 7. replicative fitness exponent
# 8. output file path

with open(argv[1]) as istream:
    ymldoc = yaml.safe_load(istream)
    ymldoc['input_files']['samp_scheme'] = argv[2]
    ymldoc['input_files']['founder_seqs'] = argv[3]
    ymldoc['input_files']['epitope_locations'] = argv[4]
    ymldoc['input_files']['conserved_sites'] = argv[5]
    ymldoc['input_files']['ref_seq'] = argv[6]
    ymldoc['parameters']['mut_rate'] = float(argv[7])
    ymldoc['parameters']['replicative_cost'] = float(argv[8])

with open(argv[9], 'w') as ostream:
    yaml.dump(ymldoc, ostream, default_flow_style=False, sort_keys=False)
