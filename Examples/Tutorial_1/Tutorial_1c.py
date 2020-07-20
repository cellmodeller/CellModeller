import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import math

def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, jitter_z=True, gamma = 20, max_planes=5)

    biophys.addPlane((0,0,0),(0,0,1),1.0) #Base plane
    biophys.addPlane((10,0,0),(-1,0,0),1.0)
    biophys.addPlane((-10,0,0),(1,0,0),1.0)
    biophys.addPlane((0,10,0),(0,-1,0),1.0)
    biophys.addPlane((0,-10,0),(0,1,0),1.0)

    # use this file for reg too
    regul = ModuleRegulator(sim, sim.moduleName)	
    # Only biophys and regulation
    sim.init(biophys, regul, None, None)
 
    # Specify the initial cell and its location in the simulation
    sim.addCell(cellType=0, pos=(0,0,0.5))

    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)
    sim.pickleSteps = 20

def init(cell):
    # Specify mean and distribution of initial cell size
    cell.targetVol = 3.0 + random.uniform(0.0,0.5)
    # Specify growth rate of cells
    cell.growthRate = 1.0
    cell.color = (1.0,1.0,1.0)

def update(cells):
    #Iterate through each cell and flag cells that reach target size for division
    for (id, cell) in cells.items():
        if cell.volume > cell.targetVol:
            cell.divideFlag = True

def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    d1.targetVol = 3.0 + random.uniform(0.0,0.5)
    d2.targetVol = 3.0 + random.uniform(0.0,0.5)

