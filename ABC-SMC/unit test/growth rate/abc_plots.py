import os
import matplotlib.pyplot as plt
import numpy as np
import pandas
import pyabc

plt.rcParams['font.size'] = 20
colors = ['b', 'g', 'm', 'c']


def abc_simulation(params):
    return None


def plot_1d(history, key, limits, exact):

    x_min = limits[key][0]
    x_max = limits[key][1]

    fig, ax = plt.subplots(figsize=(18, 15))

    for t in range(history.max_t + 1):
        df, w = history.get_distribution(m=0, t=t)
    pyabc.visualization.plot_kde_1d(
        df,
        w,
        xmin=x_min,
        xmax=x_max,
        x=key,
        xname=key,
        ax=ax,
        label="Iteration={}".format(t),
    )

    ax.axvline(exact[key], color="k", linestyle="dashed")
    ax.legend()

    ax.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.13))
    plt.savefig("images/%s_1D_KDE.png" % key, bbox_inches='tight')


def plot_matrix(history, exact):
    fig, ax = plt.subplots(3, 3, figsize=(22, 18))
    df, w = history.get_distribution()
    pyabc.visualization.plot_kde_matrix(df, w, arr_ax=ax, refval=exact, refval_color='r')
    plt.savefig("images/matrix_distributions.png", bbox_inches='tight')


def plot_sample_number(history):
    fig, ax = plt.subplots()
    pyabc.visualization.plot_sample_numbers(history, ax=ax)
    plt.savefig('images/sample_numbers.png', bbox_inches="tight")


def plot_epsilons(history):
    fig, ax = plt.subplots(figsize=(10, 7))
    pyabc.visualization.plot_epsilons(history, ax=ax)
    # pyabc.visualization.plot_effective_sample_sizes(history, ax=arr_ax[2])
    plt.savefig('images/epsilon.png', bbox_inches="tight")


def run_plotting_abc(param_config_file, exact_values, run_id):

    # generate the parameter prior distribution
    prior_distributions = {}
    # parameter value limits used in ABC computation
    limits = {}
    # param_config = open(args.param_config).read().splitlines()
    param_config = open(param_config_file).read().splitlines()

    for param_info in param_config:
        # split line in parameter file to get name, low, hi
        split_line = param_info.split(",")
        param_name = split_line[0]
        param_low = float(split_line[1])
        param_hi = float(split_line[2])
        width = abs(param_hi - param_low)
        limits[param_name] = [param_low, param_hi]

        # generate dictionary containing keywords for pyabc Distribution() class
        prior_distributions[param_name] = {"type": "uniform", "kwargs": {"loc": param_low + 1e-2, "scale": width}}

    parameter_prior = pyabc.Distribution()
    parameter_prior = parameter_prior.from_dictionary_of_dictionaries(prior_distributions)

    # create "Null" ABCSMC object that has the corresponding parameter priors used in ABC run 
    # needed to load the run history from db file
    abc = pyabc.ABCSMC(abc_simulation, parameter_prior)

    # db file created or reused during ABC computation
    db_path = ("sqlite:///" + os.getcwd() + "/spring_constants_run.db")
    print(db_path)

    # the run id corresponding to the ABC run of interest in db file

    # load history of ABC run=run_id from the db file
    history = abc.load(db_path, run_id)

    # create a matrix plot of each pairwise 2D distribution 
    plot_matrix(history, exact_values)

    # plot epsilon schedule
    plot_epsilons(history)

    # plot sample number
    plot_sample_number(history)

    # plot the 1D kde distributions for each parameter
    for key in exact_values.keys():
        plot_1d(history, key, limits, exact_values)


if __name__ == '__main__':
    param_config_file = '../../input/param_config_growth_rate.txt'
    exact_values = {"growth_rate": 0.02}
    # the run id corresponding to the ABC run of interest in db file
    run_id = 1
    run_plotting_abc(param_config_file, exact_values, run_id)

