import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
import numpy

max_cells = 2**15

#  This example shows how to load simulation data, and then apply 
#  another model to simulate more time steps.
#
#  This can be useful for example if you want to examine the properties
#  of a large colony _after_ it has reached some given size.
#

def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, \
                            max_cells=max_cells, \
                            max_contacts=32, \
                            max_sqs=128**2, \
                            jitter_z=False, \
                            gamma=10)

    biophys.addPlane((0,0,-0.5), (0,0,1), 0.2)
    biophys.addPlane((0,0,0.5), (0,0,-1), 1e-4)

    regul = ModuleRegulator(sim, sim.moduleName)	# use this file for reg too
    # Only biophys and regulation
    sim.init(biophys, regul, None, None)

    #load data from pickle file
    dataFileName = '' # Edit to include your filename here

    # Or pop up a dialog to choose a pickle
    if dataFileName=='':
        print("Please edit the model file to specify a pickle file.") 
        print("  -- No data was loaded and there are no cells in this simulation!")
    else:
        # Import the data and load into Simulator
        print("Loading data from pickle file: %s"%(dataFileName))
        import pickle
        data = pickle.load(open(dataFileName,'r'))
        sim.loadFromPickle(data)

    if sim.is_gui:
        # Add some objects to draw the models
        from CellModeller.GUI import Renderers
        therenderer = Renderers.GLBacteriumRenderer(sim)
        sim.addRenderer(therenderer)

def init(cell):
    cell.targetVol = 2.5 + random.uniform(0.0,0.5)
    cell.growthRate = 0.5

def update(cells):
    for (id, cell) in cells.items():
        cell.color = [1,0,0]
        if cell.volume > cell.targetVol:
            cell.divideFlag = True

def divide(parent, d1, d2):
    d1.targetVol = 2.5 + random.uniform(0.0,0.5)
    d2.targetVol = 2.5 + random.uniform(0.0,0.5)

