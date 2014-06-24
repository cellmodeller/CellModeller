"""This module holds the default CellModeller Parameters.

The idea is to have a single place where all default CellModeller parameters
can be seen.
The object is intended to be used as static class.

Users import this class in their model file and change the parameters
or add new ones for use in their module.

The parameters, as they stand when Cell Modeller makes the data directory, 
are recorded in comments at the start of copy of the model file. 

"""


class Parameters(object):

    """ These are the default parameters used by CellModeller.

    You may override these parameters in a user defined model.

    """

    bio_physics_dt = 0.025
    integration_dt = 0.025
    target_cell_length = 3.0 # micro meters
    cell_radius = 0.25 # micro meters
    growthRate = 0.5
    open_cl_device = 0 
    open_cl_platform = 0
