import numpy as np
import pandas as pd
from numpy.linalg import norm
from scipy.stats import norm


def fit_normal_distribution(feature):
    """
    @param feature list  feature values list
    """
    mu, std = norm.fit(feature)

    return mu

def calc_summary_statistics_for_exp(exp_data):

    """
    @param sim_data dataframe   image data for ENTIRE simulation
    # note that: you should use output of https://github.com/ingallslab/ImageProcessing as input for this function
    return: mu     float        mean of distribution
    """

    # RFP: cell type = 1, YFP: cell type = 2
    yfp_bac = exp_data.loc[exp_data['cellType'] == 2].reset_index(drop=True)
    # remove growth rate = Nan
    yfp_bac = yfp_bac.loc[(yfp_bac["growthRate"].notnull())].reset_index(drop=True)
    growth_rate = yfp_bac['growthRate'].values.tolist()

    # fit normal distribution
    mu = fit_normal_distribution(growth_rate)

    return mu


def calc_summary_statistics_for_simulation(sim_data, interval_time):
    """
    @param exp_data dataframe   image data for simulation data
    return: mu     float        mean of distribution
    """

    checked_id = []
    growth_rate = []

    for index, row in sim_data.iterrows():
        if sim_data.iloc[index]['Id'] not in checked_id:
            if sim_data.iloc[index]["YFP"] == 1:
                same_df = sim_data.loc[
                    (sim_data['Id'] == sim_data.iloc[index]["Id"])].reset_index(drop=True)
                same_length = same_df['AreaShape_MajorAxisLength'].values.tolist()
                if len(same_length) > 1 and same_length[-1] > 0 and same_length[0] > 0:
                    t = len(same_length) * interval_time
                    elongation_rate = round((np.log(same_length[-1]) - np.log(same_length[0])) / t, 3)
                    growth_rate.append(elongation_rate)
            checked_id.append(sim_data.iloc[index]["Id"])

    # fit normal distribution
    mu = fit_normal_distribution(growth_rate)

    return mu


def get_simulation_summary_statistics(sim_data, interval_time):
    """
    Main function for getting feature data and packaging each feature data
    object into a dictionary for later use.

    Inputs:
        sim_data (pandas dataframe) -> image data for ENTIRE simulation
        image_count (int) -> image number
    Outputs:
        feature_dict (dict) --> Dictionary containing keyed entries for all
                                features and their data objects
    """
    features = calc_summary_statistics_for_simulation(sim_data, interval_time)
    # build feature data structure
    feature_dict = {"mean growth rate": features}

    return feature_dict


def get_exp_summary_statistics(exp_data):
    """
    @param exp_data   dataframe   image data for observation data
    return feature_dict  dict     Dictionary containing keyed entries for all features and their data objects
    """

    features = calc_summary_statistics_for_exp(exp_data)
    # build feature data structure
    feature_dict = {"mean growth rate": features}

    return feature_dict
