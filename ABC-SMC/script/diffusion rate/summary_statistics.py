import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix
from sklearn.linear_model import LinearRegression


def fit_linear_regression(all_yfp_in_current_timestep_distance, delta_l_list):
    # create object for the class
    linear_regressor = LinearRegression()

    linear_regressor.fit(np.array(all_yfp_in_current_timestep_distance).reshape(-1, 1),
                         np.array(delta_l_list).reshape(-1, 1))  # perform linear regression
    intercept = linear_regressor.intercept_[0]
    slope = linear_regressor.coef_[0][0]
    line_space_for_linear_regression = np.linspace(min(all_yfp_in_current_timestep_distance),
                                                   max(all_yfp_in_current_timestep_distance), 500)

    return line_space_for_linear_regression, slope, intercept


def calc_summary_statistics_for_simulation(sim_data):
    """
    @param sim_data dataframe   image data for ENTIRE simulation
    return: time_step_list, line_space_for_linear_regression_list, slope_list, intercept_list    list
    """

    time_step_list = []
    line_space_for_linear_regression_list = []
    slope_list = []
    intercept_list = []

    # loop over images in Sim data
    image_number_list = sim_data["ImageNumber"].unique()  # [1,2,3,4,5,....,30]

    for image_number in image_number_list:

        yfp_delta_length = []
        yfp_distance_from_rfp_in_current_timestep = []

        # obtain pandas dataframe containing data corresponding to
        current_df = sim_data.loc[sim_data["ImageNumber"] == image_number].reset_index(drop=True)

        # obtain pertinent position data from dataframe
        current_rfp_bac = current_df.loc[current_df["RFP"] == 1].reset_index(drop=True)
        current_yfp_bac = current_df.loc[current_df["YFP"] == 1].reset_index(drop=True)

        if current_rfp_bac.shape[0] != 0 and current_yfp_bac.shape[0] != 0:
            distance_yfp_rfp_df = pd.DataFrame(
                distance_matrix(current_yfp_bac[["AreaShape_Center_X", "AreaShape_Center_Y"]].values,
                                current_rfp_bac[["AreaShape_Center_X", "AreaShape_Center_Y"]].values),
                index=current_yfp_bac.index, columns=current_rfp_bac.index)

            for index, row in current_yfp_bac.iterrows():
                # 1/(1+di)^2
                distance_val = (1 / (1 + distance_yfp_rfp_df.iloc[index]) ** 2)
                sum_of_distance_from_rfp_bacteria = distance_val.sum()
                # calculate delta Length
                next_timestep = image_number + 1
                # same bacteria
                same_bac_in_next_timestep = sim_data.loc[(sim_data["ImageNumber"] == next_timestep) &
                                                         (sim_data['Id'] == current_yfp_bac.iloc[index]['Id'])
                                                         ].reset_index(drop=True)
                # it means bacterium life history has been continued
                if same_bac_in_next_timestep.shape[0] != 0:
                    delta_length = same_bac_in_next_timestep['AreaShape_MajorAxisLength'].values.tolist()[0] - \
                                   current_yfp_bac.iloc[index]["AreaShape_MajorAxisLength"]

                    yfp_distance_from_rfp_in_current_timestep.append(sum_of_distance_from_rfp_bacteria)
                    yfp_delta_length.append(delta_length)

        # fit linear regression
        if len(yfp_distance_from_rfp_in_current_timestep) > 0:
            line_space_for_linear_regression, slope, intercept = \
                fit_linear_regression(yfp_distance_from_rfp_in_current_timestep, yfp_delta_length)

            time_step_list.append(image_number)
            line_space_for_linear_regression_list.append(line_space_for_linear_regression)
            slope_list.append(slope)
            intercept_list.append(intercept)

    return time_step_list, line_space_for_linear_regression_list, slope_list, intercept_list


def calc_summary_statistics_for_exp(exp_data):
    """
    @param exp_data dataframe   image data for observation data
    # note that: you should use output of https://github.com/ingallslab/ImageProcessing as input for this function
    return: line_space_for_linear_regression, slope, intercept    list
    """

    time_step_list = []
    line_space_for_linear_regression_list = []
    slope_list = []
    intercept_list = []

    # loop over images in Sim data
    image_number_list = exp_data["stepNum"].unique()

    for image_number in image_number_list:

        yfp_delta_length_list = []
        yfp_distance_from_rfp_in_current_timestep = []

        # obtain pandas dataframe containing data corresponding to
        current_df = exp_data.loc[exp_data["stepNum"] == image_number].reset_index(drop=True)

        # obtain pertinent position data from dataframe
        # RFP: cell type = 1, YFP: cell type = 2
        current_rfp_bac = current_df.loc[current_df["cellType"] == 1].reset_index(drop=True)
        current_yfp_bac = current_df.loc[current_df["cellType"] == 2].reset_index(drop=True)

        if current_rfp_bac.shape[0] != 0 and current_yfp_bac.shape[0] != 0:
            distance_yfp_rfp_df = pd.DataFrame(
                distance_matrix(current_yfp_bac[["X", "Y"]].values, current_rfp_bac[["X", "Y"]].values),
                index=current_yfp_bac.index, columns=current_rfp_bac.index)

            for index, row in current_yfp_bac.iterrows():
                # 1/(1+di)^2
                distance_val = (1 / (1 + distance_yfp_rfp_df.iloc[index]) ** 2)
                sum_of_distance_from_rfp_bacteria = distance_val.sum()
                # calculate delta L
                next_timestep = image_number + 1
                # find same bacteria
                same_bac_in_next_timestep = exp_data.loc[(exp_data["stepNum"] == next_timestep) &
                                                         (exp_data['id'] == current_yfp_bac.iloc[index]["id"])
                                                         ].reset_index(drop=True)
                if same_bac_in_next_timestep.shape[0] != 0:
                    delta_length = same_bac_in_next_timestep['length'].values.tolist()[0] - \
                                   current_yfp_bac.iloc[index]["length"]

                    yfp_distance_from_rfp_in_current_timestep.append(sum_of_distance_from_rfp_bacteria)
                    yfp_delta_length_list.append(delta_length)

        # fit linear regression
        if len(yfp_distance_from_rfp_in_current_timestep) > 0:
            line_space_for_linear_regression, slope, intercept = \
                fit_linear_regression(yfp_distance_from_rfp_in_current_timestep, yfp_delta_length_list)

            time_step_list.append(image_number)
            line_space_for_linear_regression_list.append(line_space_for_linear_regression)
            slope_list.append(slope)
            intercept_list.append(intercept)

    return time_step_list, line_space_for_linear_regression_list, slope_list, intercept_list


def get_exp_summary_statistics(exp_data):
    """
    @param exp_data   dataframe   image data for observation data
    return feature_dict  dict     Dictionary containing keyed entries for all features and their data objects
    """

    time_step, line_space, slope, intercept = calc_summary_statistics_for_exp(exp_data)

    list_to_dct = {'time step list': time_step, 'line space list': line_space, 'slope list': slope,
                   'intercept list': intercept}

    return list_to_dct


def get_simulation_summary_statistics(sim_data):

    """
    @param sim_data         dataframe       image data for ENTIRE simulation
    Outputs: feature_dict   dict            Dictionary containing keyed entries for all features and their data objects
    """

    time_step, line_space, slope, intercept = calc_summary_statistics_for_simulation(sim_data)

    list_to_dct = {'time step list': time_step, 'line space list': line_space, 'slope list': slope,
                   'intercept list': intercept}

    return list_to_dct
