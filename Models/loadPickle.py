import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy

max_cells = 2**15

def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, max_cells=max_cells, max_contacts=32, max_sqs=64**2, jitter_z=False, reg_param=2)
    
    regul = ModuleRegulator(sim, __file__)	# use this file for reg too

    sim.init(biophys, regul, None, None)
    
    # load from pickle file
    import pickle
    (cs, lineage) = pickle.load(open('data/biophysCL-01-03-13-04-12/step-02450.pickle','r'))
    print "Loading %i cells..."%len(cs)
    sim.setCellStates(cs)
    biophys.load_from_cellstates(cs)

    #for (i,cell) in sim.cellStates.items():
    #    cell.color = [0.1, 0.1+1000.0*cell.species[3], 
    #    0.1+200.0*cell.species[4]]
#    for (i,cell) in sim.cellStates.items():
#        cell.color = [0.1, cell.species[3]*0.1, cell.species[4]*0.1]

    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)

def init(cell):
    pass
    
def numSignals():
    return 0  

def numSpecies():
    return 0

def signalRates(cell, speciesLevels, signalLevels):
    return []

def speciesRates(cell, speciesLevels, signalLevels):
    return []
    
def update(cells):
    pass
