import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import math

rodlen = 4.0

max_cells = 400000

#cell_colors = {0:[0.0, 1.0, 0.0],
#               1:[0.0, 0.0, 1.0],
#               2:[1.0, 0.0, 0.0],
#               3:[0.0, 1.0, 1.0]}
cell_colors = numpy.random.uniform(0,1,(9,3))


def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, max_substeps=8, max_cells=max_cells, max_contacts=32, max_sqs=192**2, jitter_z=False, reg_param=0.04, gamma=30)
    #biophys.addPlane((0,0,-0.5), (0,0,1), 1.0)
    #biophys.addPlane((0,0,0.5), (0,0,-1), math.sqrt(7.5e-4))


    regul = ModuleRegulator(sim, __file__)	# use this file for reg too
    # Only biophys and regulation
    sim.init(biophys, regul, None, None)

    sim.addCell(cellType=0, len=rodlen, pos=(0,0,0))

    #sim.addCell(cellType=0, pos=(0,-10.0,0))
    #sim.addCell(cellType=1, pos=(0,10.0,0))

    #sim.addCell(cellType=0, pos=(16,16,0))
    #sim.addCell(cellType=1, pos=(0,16,0))
    #sim.addCell(cellType=2, pos=(-16,16,0))
    #sim.addCell(cellType=3, pos=(16,0,0))
    #sim.addCell(cellType=4, pos=(0,0,0))
    #sim.addCell(cellType=5, pos=(-16,0,0))
    #sim.addCell(cellType=6, pos=(16,-16,0))
    #sim.addCell(cellType=7, pos=(0,-16,0))
    #sim.addCell(cellType=8, pos=(-16,-16,0))


    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)
    sim.savePickle = True
    sim.pickleSteps = 20

def init(cell):
    cell.targetVol = rodlen + random.uniform(0.0,rodlen*0.125)
    cell.growthRate = 1.0 

def numSignals():
    return 0

def numSpecies():
    return 0

def update(cells):
    for (id, cell) in cells.iteritems():
        cell.color = cell_colors[cell.cellType]
        if cell.volume > cell.targetVol:
            cell.asymm = [1,1]
            cell.divideFlag = True

def divide(parent, d1, d2):
    d1.targetVol = rodlen + random.uniform(0.0,rodlen*0.125)
    d2.targetVol = rodlen + random.uniform(0.0,rodlen*0.125)

