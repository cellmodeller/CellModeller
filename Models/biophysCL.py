import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import math

max_cells = 400000

def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, max_substeps=8, max_cells=max_cells, max_contacts=32, max_sqs=192**2, jitter_z=False, reg_param=2, gamma=10)
    #biophys.addPlane((0,0,-0.5), (0,0,1), 1.0)
    #biophys.addPlane((0,0,0.5), (0,0,-1), math.sqrt(7.5e-4))


    regul = ModuleRegulator(sim, __file__)	# use this file for reg too
    # Only biophys and regulation
    sim.init(biophys, regul, None, None)

    #ct = 0
    #for x in range(-4,5):
    #    for y in range(-4,5):
    sim.addCell(cellType=0, pos=(0,0,0)) #x*40,y*40,0))
    #        ct += 1
    #sim.addCell(cellType=1, pos=(20,0,0))


    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)
    sim.pickleSteps = 10

def init(cell):
    cell.targetVol = 2.5 + random.uniform(0.0,0.5)
    cell.growthRate = 2.0

def numSignals():
    return 0

def numSpecies():
    return 0

def update(cells):
    for (id, cell) in cells.iteritems():
        cell.color = [cell.cellType*0.6+0.1, 1.0-cell.cellType*0.6, 0.3]
        #max(cell.startVol*2.0,0.0): #cell.startvol*2: #random.uniform(1.75,2.0):
        if cell.volume > cell.targetVol:
            a = 1#random.uniform(0.95,1.05)
            cell.asymm = [a,1]
            cell.divideFlag = True

def divide(parent, d1, d2):
    d1.targetVol = 2.5 + random.uniform(0.0,0.5)
    d2.targetVol = 2.5 + random.uniform(0.0,0.5)

