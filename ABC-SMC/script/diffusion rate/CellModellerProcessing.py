import math
import numpy as np
import pickle
import CellModeller
import pandas as pd
import glob


def extract_features(current_time_step, current_bacteria_info, previous_bacteria, dataframe, bacteria_id_label):
    # bacteria information
    cs = current_bacteria_info['cellStates']

    # Cellmodeller pickle file structure:
    # ['cellStates', 'stepNum', 'lineage', 'moduleStr', 'moduleName']
    # important keys: 1: 'cellStates' (type: dictionary)  2: 'lineage' (dictionary)
    # 'cellStates' dictionary keys: bacteria id
    # 'lineage' dictionary: daughter id: parent id

    for CellID, it in enumerate(cs.keys()):
        if cs[it].id in bacteria_id_label.keys():  # it means: Life has continued for bacterium

            # last occurrence of element in list
            last_occurrence_index_in_list = max(idx for idx, val in enumerate(dataframe['Id']) if val == cs[it].id)

            # cell age
            dataframe['CellAge'].append(dataframe['CellAge'][last_occurrence_index_in_list])

            # parent information
            dataframe['TrackObjects_ParentImageNumber_50'].append(
                dataframe['ImageNumber'][last_occurrence_index_in_list])
            dataframe['TrackObjects_ParentObjectNumber_50'].append(dataframe['ObjectNumber']
                                                                   [last_occurrence_index_in_list])
            # assign label
            dataframe["TrackObjects_Label_50"].append(bacteria_id_label[cs[it].id])

        else:  # it means: A bacterium has been born or a cell division has taken place
            if current_time_step == 1:  # birth
                if len(bacteria_id_label.values()) == 0:
                    this_bacterium_label = 1
                else:
                    this_bacterium_label = max(bacteria_id_label.values()) + 1
                # cell age
                dataframe['CellAge'].append(0)
                # parent information
                dataframe['TrackObjects_ParentImageNumber_50'].append(0)
                dataframe['TrackObjects_ParentObjectNumber_50'].append(0)

            else:
                # find parent id
                if current_bacteria_info['lineage'][cs[it].id]:  # parent bacteria has been found
                    parent_id = current_bacteria_info['lineage'][cs[it].id]
                    # assign label
                    if parent_id in bacteria_id_label.keys():
                        this_bacterium_label = bacteria_id_label[parent_id]
                        # last occurrence of element in list
                        last_occurrence_index_in_list = max(
                            idx for idx, val in enumerate(dataframe['Id']) if val == parent_id)
                    else:
                        '''
                            Find the nearest bacterium from the previous timestep that does not exist in the current
                             timestep (since the founded bacterium is the grandmother of the daughter)                    
                        '''
                        # find Bacteria whose life has ended.
                        life_ended_bacteria_in_previous_time_step = [bacterium_id for bacterium_id in
                                                                     previous_bacteria['lineage'].keys() if
                                                                     bacterium_id not in current_bacteria_info[
                                                                         'lineage'].values()]
                        previous_bacteria_info = previous_bacteria['cellStates']
                        # calculate distance
                        distance = []
                        for life_ended_bacteria_id in life_ended_bacteria_in_previous_time_step:
                            this_bacterium_center = np.array((cs[it].pos[0], cs[it].pos[1]))
                            life_ended_bacterium_in_previous_time_step = np.array(
                                (previous_bacteria_info[life_ended_bacteria_id].pos[0],
                                 previous_bacteria_info[life_ended_bacteria_id].pos[1]))
                            # find distance
                            distance_value = np.linalg.norm(
                                this_bacterium_center - life_ended_bacterium_in_previous_time_step)
                            distance.append(distance_value)
                        # find the nearest grandmother bacterium
                        grandmother_id = sorted(zip(distance, life_ended_bacteria_in_previous_time_step))[0][1]
                        # assign label
                        this_bacterium_label = bacteria_id_label[grandmother_id]
                        last_occurrence_index_in_list = max(
                            idx for idx, val in enumerate(dataframe['Id']) if val == grandmother_id)

                else:
                    '''
                        Find the nearest bacterium from the previous timestep that does not exist in the current timestep
                         (since the founded bacterium is the grandmother of the daughter)                    
                    '''
                    # find Bacteria whose life has ended.
                    life_ended_bacteria_in_previous_time_step = [bacterium_id for bacterium_id in
                                                                 previous_bacteria['lineage'].keys() if
                                                                 bacterium_id not in current_bacteria_info[
                                                                     'lineage'].values()]
                    previous_bacteria_info = previous_bacteria['cellStates']
                    # calculate distance
                    distance = []
                    for life_ended_bacteria_id in life_ended_bacteria_in_previous_time_step:
                        this_bacterium_center = np.array((cs[it].pos[0], cs[it].pos[1]))
                        life_ended_bacterium_in_previous_time_step = np.array(
                            (previous_bacteria_info[life_ended_bacteria_id].pos[0],
                             previous_bacteria_info[life_ended_bacteria_id].pos[1]))
                        # find distance
                        distance_value = np.linalg.norm(
                            this_bacterium_center - life_ended_bacterium_in_previous_time_step)
                        distance.append(distance_value)
                    # find the nearest grandmother bacterium
                    grandmother_id = sorted(zip(distance, life_ended_bacteria_in_previous_time_step))[0][1]
                    # assign label
                    this_bacterium_label = bacteria_id_label[grandmother_id]
                    last_occurrence_index_in_list = max(
                        idx for idx, val in enumerate(dataframe['Id']) if val == grandmother_id)

                # cell age
                dataframe['CellAge'].append(dataframe['CellAge'][last_occurrence_index_in_list] + 1)

                # parent information
                dataframe['TrackObjects_ParentImageNumber_50'].append(
                    dataframe['ImageNumber'][last_occurrence_index_in_list])
                dataframe['TrackObjects_ParentObjectNumber_50'].append(dataframe['ObjectNumber']
                                                                       [last_occurrence_index_in_list])

            # assign label
            dataframe["TrackObjects_Label_50"].append(this_bacterium_label)
            # add new key: value to dictionary
            bacteria_id_label[cs[it].id] = this_bacterium_label

        dataframe['Id'].append(cs[it].id)
        dataframe['ImageName'].append(current_bacteria_info['stepNum'])
        dataframe['ImageNumber'].append(current_time_step)
        dataframe['ObjectNumber'].append(CellID + 1)

        # cell Type
        # In CellModeller CellTypes are: 0,1,2,3,...
        for index in range(len(dataframe['Type'])):
            if cs[it].cellType == index:
                dataframe['Type'][index].append(1)
            else:
                dataframe['Type'][index].append(0)

        dataframe['AreaShape_Center_X'].append(cs[it].pos[0])
        dataframe['AreaShape_Center_Y'].append(cs[it].pos[1])
        dataframe['AreaShape_MajorAxisLength'].append(cs[it].length)
        dataframe['AreaShape_MinorAxisLength'].append(cs[it].radius)
        angle = math.atan2(cs[it].dir[1], cs[it].dir[0])
        dataframe['AreaShape_Orientation'].append(angle)
        dataframe['Node_x1_x'].append(cs[it].ends[0][0])
        dataframe['Node_x1_y'].append(cs[it].ends[0][1])
        dataframe['Node_x2_x'].append(cs[it].ends[1][0])
        dataframe['Node_x2_y'].append(cs[it].ends[1][1])
        # Surface area of a capsule:
        # S = 2πr(2r + a)
        dataframe['AreaShape_Area'].append(2 * math.pi * cs[it].radius * (cs[it].length + 2 * cs[it].radius))

    return dataframe, bacteria_id_label


def starting_process(input_directory, cell_types):
    """

    @param input_directory  str  directory of pickle files:
    @param cell_types       list  bacteria cell types list.
    The list of cell types must match the list defined in the CellModeller simulation file.
    (e.g: In the case where cell type 1 shows RFP, and cell type 2 shows YFP, the cell type list should be:
     ['RFP', 'YFP'])
    @return: df dataframe

    """

    # firstly I create a dictionary and append extracted features to corresponding key list

    dataframe = {'Id': [], 'ImageNumber': [], 'ObjectNumber': [], 'Type': [], 'AreaShape_Area': [],
                 'AreaShape_Center_X': [], 'AreaShape_Center_Y': [], 'AreaShape_MajorAxisLength': [],
                 'AreaShape_MinorAxisLength': [], 'AreaShape_Orientation': [], 'Node_x1_x': [], 'Node_x1_y': [],
                 'Node_x2_x': [], 'Node_x2_y': [], 'CellAge': [], 'TrackObjects_ParentImageNumber_50': [],
                 'TrackObjects_ParentObjectNumber_50': [], 'validID': [], 'ImageName': [], 'TrackObjects_Label_50': []}

    # keys: bacteria id
    # values: assigned bacteria labels
    bacteria_id_label = {}

    for cell_type in cell_types:
        dataframe['Type'].append([])

    # read pickle files
    path = input_directory + "/*.pickle"
    filename_list = [filename for filename in sorted(glob.glob(path))]

    for cnt, filename in enumerate(filename_list):

        # read current pickle file
        current_bacteria_info = pickle.load(open(filename, 'rb'))
        time_step = cnt + 1

        # read previous pickle file
        if cnt > 0:
            previous_bacteria = pickle.load(open(filename_list[cnt - 1], 'rb'))
        else:
            previous_bacteria = ''

        # extract features
        dataframe, bacteria_id_label = extract_features(time_step, current_bacteria_info, previous_bacteria,
                                                        dataframe, bacteria_id_label)

    # create data frame
    df = pd.DataFrame({'Id': dataframe['Id'], 'ImageName': dataframe['ImageName'],
                       'ImageNumber': dataframe['ImageNumber'], 'ObjectNumber': dataframe['ObjectNumber'],
                       'AreaShape_Area': dataframe['AreaShape_Area'],
                       'AreaShape_Center_X': dataframe['AreaShape_Center_X'],
                       'AreaShape_Center_Y': dataframe['AreaShape_Center_Y'],
                       'AreaShape_MajorAxisLength': dataframe['AreaShape_MajorAxisLength'],
                       'AreaShape_MinorAxisLength': dataframe['AreaShape_MinorAxisLength'],
                       'AreaShape_Orientation': dataframe['AreaShape_Orientation'],
                       "TrackObjects_Label_50": dataframe["TrackObjects_Label_50"],
                       'TrackObjects_ParentImageNumber_50': dataframe['TrackObjects_ParentImageNumber_50'],
                       'TrackObjects_ParentObjectNumber_50': dataframe['TrackObjects_ParentObjectNumber_50'],
                       'CellAge_Generation': dataframe['CellAge'], 'Node_x1_x': dataframe['Node_x1_x'],
                       'Node_x1_y': dataframe['Node_x1_y'], 'Node_x2_x': dataframe['Node_x2_x'],
                       'Node_x2_y': dataframe['Node_x2_y']})

    # Using DataFrame.insert() to add a column
    # add new columns after third column
    column_num = 3
    for cnt, CellType in enumerate(cell_types):
        df.insert(column_num, CellType, dataframe['Type'][cnt], True)
        column_num += 1

    return df
