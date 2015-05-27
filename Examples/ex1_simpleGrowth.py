import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import math

cell_cols = {   0:[0.7,0.2,0], \
                1:[0,0.2,0.7], \
                2:[0.2,0.7,0], \
                3:[0.7,0.2,0.7], \
                4:[0.2,0.7,0.7]}

def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, jitter_z=False)

    # use this file for reg too
    regul = ModuleRegulator(sim, sim.moduleName)	
    # Only biophys and regulation
    sim.init(biophys, regul, None, None)
 

    # Specify the initial cell and its location in the simulation
    for i in range(5):
        px = random.uniform(-10.0,10.0)
        py = random.uniform(-10.0,10.0)
        dx = random.uniform(0,1.0)
        dy = math.sqrt(1-dx*dx)
        sim.addCell(cellType=i, pos=(px,py,0), dir=(dx,dy,0)) 

    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)
    sim.pickleSteps = 10

def init(cell):
    # Specify mean and distribution of initial cell size
    cell.targetVol = 3.5 + random.uniform(0.0,0.5)
    # Specify growth rate of cells
    cell.growthRate = 1.0

def update(cells):
    #Iterate through each cell and flag cells that reach target size for division
    for (id, cell) in cells.iteritems():
        cell.color = cell_cols[cell.cellType]
        if cell.volume > cell.targetVol:
            cell.divideFlag = True

def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    d1.targetVol = 3.5 + random.uniform(0.0,0.5)
    d2.targetVol = 3.5 + random.uniform(0.0,0.5)

