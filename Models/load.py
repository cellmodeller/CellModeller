import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy

max_cells = 400000

def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, max_cells=max_cells, max_contacts=32, max_sqs=128**2, jitter_z=False, reg_param=2, gamma=10)
    biophys.addPlane((0,0,-0.5), (0,0,1), 0.2)
    biophys.addPlane((0,0,0.5), (0,0,-1), 1e-4)


    regul = ModuleRegulator(sim, __file__)	# use this file for reg too
    # Only biophys and regulation
    sim.init(biophys, regul, None, None)

    #sim.addCell(cellType=0, pos=(0,0,0))

    #load data from pickle file
    import cPickle
    data = cPickle.load(open('/Users/tjr34/LazyLowRe/CellModeller/build/lib/data/lobes-16-06-30-08-12/step-00530.pickle','r'))
    cellStates = data[0]
    print "Loading %i cells"%len(cellStates)
    sim.cellStates = cellStates
    regul.cellStates = cellStates
    biophys.load_from_cellstates(cellStates)
    idx_map = {}
    for id,state in sim.cellStates.iteritems():
        idx_map[state.id] = state.idx
    sim.idToIdx = idx_map

#    parents = data[1]
#    for id,state in sim.cellStates.iteritems():
#        p = lineage(parents, [1,2], id)
#        if p == 1: state.color = [0, 1, 0]
#        elif p == 2: state.color = [0, 0, 1]

    #import numpy.random
    #cellcols = numpy.random.uniform(0,1,81*3)
    #cellcols.shape = (81,3)
    #for id,state in sim.cellStates.iteritems():
    #    state.color = list(cellcols[state.cellType,:])

    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)


def lineage(parents, founders, id):
    while id not in founders:
        id = parents[id]
    return id


def init(cell):
    cell.targetVol = 2.5 + random.uniform(0.0,0.5)
    cell.growthRate = 0.5

def numSignals():
    return 0

def numSpecies():
    return 0

def update(cells):
    pass

def divide(parent, d1, d2):
    d1.targetVol = 2.5 + random.uniform(0.0,0.5)
    d2.targetVol = 2.5 + random.uniform(0.0,0.5)

