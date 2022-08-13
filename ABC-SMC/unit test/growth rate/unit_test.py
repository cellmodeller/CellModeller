import os
import sys
import pandas as pd
sys.path.append('../../script/growth rate')
from Running_ABC_SMC import running_abc

# Parameter list with containing lines with format
# '<parameter name>,<low value>,<high value>
param_config = '../../input/param_config_growth_rate.txt'
# dataset used to compare ABC simulations against
exp_data = "../../input/Ati diffusible-22-07-04-13-09.csv"

output_directory = os.getcwd()
print(output_directory)

# configure ABC specific arguments
number_of_particles = 8
max_populations = 2

# read txt file
param_config = open(param_config).read().splitlines()

# read experimental data
exp_data = pd.read_csv(exp_data)

running_abc(exp_data, param_config, output_directory, number_of_particles, max_populations)
