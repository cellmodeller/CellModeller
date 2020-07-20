import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy as np
import math

max_cells = 50000


def setup(sim):
    # Set biophysics, signalling, and regulation models
    
    #max_sqr
    #biophys = CLBacterium(sim, max_cells=max_cells, jitter_z=False, max_sqs=256**2)
    biophys = CLBacterium(sim, max_cells=max_cells, jitter_z=False, max_spheres=1, gamma=100.)

    # use this file for reg too
    regul = ModuleRegulator(sim)
    # Only biophys and regulation
    sim.init(biophys, regul, None, None)

    # Specify the initial cell and its location in the simulation
    sim.addCell(cellType=0, pos=(0,0,0))

    biophys.addSphere((0,0,0), 10, 1.0, -1)


    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)

    sim.pickleSteps = 5 

def init(cell):
    # Specify mean and distribution of initial cell size
    cell.targetVol = 3.0 + random.uniform(0.0,0.5)
    # Specify growth rate of cells
    #cell.growthRate = 0.6
    cell.growthRate = 1.0
    cell.color = [0, 1, 0]

def update(cells):
    #Iterate through each cell and flag cells that reach target size for division
    for (id, cell) in cells.items():
        if cell.volume > cell.targetVol:
            cell.divideFlag = True

def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    d1.targetVol = 3.0 + random.uniform(0.0,0.5)
    d2.targetVol = 3.0 + random.uniform(0.0,0.5)

