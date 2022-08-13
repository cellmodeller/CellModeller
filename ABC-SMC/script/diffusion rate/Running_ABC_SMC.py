# Approximate Bayesian Computation with pyabc
import pyabc
import sys

sys.path.append('../../../Scripts/')
from CellModeller.Simulator import Simulator
import os
import subprocess
import multiprocessing
import numpy as np
import pandas
import random
import string
from CreateSimulationScript import simulation_script
from CellModellerProcessing import starting_process
from summary_statistics import get_exp_summary_statistics, get_simulation_summary_statistics


def get_features(bacteria_simulation_df):
    cell_modeller_features = get_simulation_summary_statistics(bacteria_simulation_df)

    return cell_modeller_features


# distance_function
def distance_calculation(cell_modeller_data_summary_statistics, exp_summary_statistics):
    distance_sim_exp = []

    for i, time_step in enumerate(cell_modeller_data_summary_statistics['time step list']):
        try:
            exp_same_time_step_index = exp_summary_statistics['time step list'].index(time_step)
        except:
            # not found
            exp_same_time_step_index = -1

        if exp_same_time_step_index != -1:
            lower_bound = max(cell_modeller_data_summary_statistics['line space list'][i][0],
                              exp_summary_statistics['line space list'][exp_same_time_step_index][0])
            upper_bound = min(cell_modeller_data_summary_statistics['line space list'][i][-1],
                              exp_summary_statistics['line space list'][exp_same_time_step_index][-1])
            x_space = np.linspace(lower_bound, upper_bound, 500)

            y_cell_modeller = cell_modeller_data_summary_statistics['slope list'][i] * x_space + \
                              cell_modeller_data_summary_statistics['intercept list'][i]
            y_exp = exp_summary_statistics['slope list'][exp_same_time_step_index] * x_space + \
                    exp_summary_statistics['intercept list'][exp_same_time_step_index]

            if exp_summary_statistics['slope list'][exp_same_time_step_index] > \
                    cell_modeller_data_summary_statistics['slope list'][i]:
                subtract = [a_i - b_i for a_i, b_i in zip(y_exp, y_cell_modeller)]
            else:
                subtract = [a_i - b_i for a_i, b_i in zip(y_cell_modeller, y_exp)]

            # distance of two lines
            area = np.trapz(y=subtract, x=x_space)
            distance_sim_exp.append(area)

    mean_distance = np.mean(distance_sim_exp)

    return mean_distance


def simulate(modfilename, platform, device, output_name, steps=50):

    total_simulation_time = 92.8
    dt = 0.016

    (path, name) = os.path.split(modfilename)
    modname = str(name).split('.')[0]
    sys.path.append(path)
    sim = Simulator(modname, 0.016, outputDirName=output_name, clPlatformNum=platform, clDeviceNum=device,
                    saveOutput=True)
    while sim.stepNum * dt <= total_simulation_time:
        sim.step()


def abc_simulation(params):
    # return loaded csv file of sim data
    # Execute CellModeller simulation with given parameters
    script_name = script_name = ''.join(
        random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(6))
    print(script_name)

    parameter = ','.join([str(params[k]) for k in params])
    simulation_script(script_name, str(parameter))

    script_path = "scripts/" + script_name + ".py"
    # sys.stdout = open(os.devnull, 'w')
    simulate(script_path, 0, 0, script_name)
    print("Finished Simulation")

    # converting CellModeller Output Structure to CP Output Structure
    # Getting the location of simulations.
    pickle_directory = "data/" + script_name
    cell_types = ['RFP', 'YFP']

    # Start Processing
    cell_modeller_df = starting_process(pickle_directory, cell_types)
    print("start post-processing")
    return get_features(cell_modeller_df)


def running_abc(exp_data, param_config, output_directory, n_particles, max_populations, epsilon):
    # create directories
    simulation_result_path = "temp/data"
    scripts_path = "temp/scripts"
    if not os.path.exists(simulation_result_path):
        os.makedirs(simulation_result_path)
    if not os.path.exists(scripts_path):
        os.makedirs(scripts_path)
    # change directory
    os.chdir("temp")

    print("Running ABC test")
    # parameter set for AstroABC package
    prior_distributions = {}

    for parameter_info in param_config:
        # split line in parameter file to get name, low, hi
        split_line = parameter_info.split(",")
        param_name = split_line[0]
        param_low = float(split_line[1])
        param_hi = float(split_line[2])
        width = abs(param_hi - param_low)

        # generate dictionary containing keywords for pyabc Distribution() class
        prior_distributions[param_name] = {"type": "uniform", "kwargs": {"loc": param_low + 1e-2, "scale": width}}

    parameter_prior = pyabc.Distribution()
    parameter_prior = parameter_prior.from_dictionary_of_dictionaries(prior_distributions)

    # extract features from experimental data
    exp_features = get_simulation_summary_statistics(exp_data)

    # running ABC
    abc = pyabc.ABCSMC(models=abc_simulation, parameter_priors=parameter_prior, distance_function=distance_calculation,
                       population_size=int(n_particles))

    db_path = ("sqlite:///" + output_directory + "/spring_constants_run.db")

    abc.new(db_path, exp_features)

    abc.run(max_nr_populations=int(max_populations), minimum_epsilon=epsilon)

    print('Finish')