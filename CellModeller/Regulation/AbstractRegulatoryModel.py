""" This file documents the form of the user-defined python module that
constitutes a genetic regulatory model.
"""

import numpy

# Required function definitions
#
def numSpecies():
""" Tell the simulator how many cell-autonomous species there are.
"""
    raise NotImplementedError()

def numSignals():
""" Tell the simulator how many non cell-autonomous signals there are.
"""
    raise NotImplementedError()

def initCell(cellState):
""" Create/set any variables in CellState for initial conditions. This is called by the
simulator to initialise cells created before the simulation is run. It is
not called when cells are created by division. See divide() below.
"""
    raise NotImplementedError()

def speciesRates(cellState, speciesLevels, signalLevels):
""" Return numpy array with d/dt of species with given species and signal levels
"""
    raise NotImplementedError()

def signalRates(cellState, speciesLevels, signalLevels):
""" Return numpy array with d/dt of signals with given species and signal levels
"""
    raise NotImplementedError()

# Optional function definitions:
#
def divide(parent, daughter1, daughter2):
""" Implement specific partitioning model of parent state variables into daughters,
this function is called by the simulator on division of a cell.

daughter1, daughter2 will have been set with the same content as the parent.
"""
    raise NotImplementedError()


