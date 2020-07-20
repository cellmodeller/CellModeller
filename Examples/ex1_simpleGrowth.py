import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
import numpy
import math

cell_cols = {0:[0,1.0,0], 1:[1.0,0,0]}
outfile = 'all.csv'

def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, jitter_z=False, max_cells=100000, gamma=100.)

    # use this file for reg too
    regul = ModuleRegulator(sim, sim.moduleName)	
    # Only biophys and regulation
    sim.init(biophys, regul, None, None)
 
    # Specify the initial cell and its location in the simulation
    sim.addCell(cellType=0, pos=(0,0,0), dir=(1,0,0))

    # Add some objects to draw the models
    #if sim.is_gui:
    from CellModeller.GUI import Renderers
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)
    #else:
    #    print("Running in batch mode: no display will be output")

    sim.pickleSteps = 100
    sim.saveOutput = True

def init(cell):
    # Specify mean and distribution of initial cell size
    cell.targetVol = 3.5 + random.uniform(0.0,0.5)
    # Specify growth rate of cells
    cell.growthRate = 1.0
    cell.color = cell_cols[cell.cellType]

def update(cells):
    #Iterate through each cell and flag cells that reach target size for division
    for (id, cell) in cells.items():
        if cell.volume > cell.targetVol:
            cell.divideFlag = True

        gr = cell.strainRate/0.05
        cgr = gr - 0.5
        # Return value is tuple of colors, (fill, stroke)
        #if cgr>0:
        #    cell.color = [1, 1-cgr*2, 1-cgr*2]
        #else:
        #    cell.color = [1+cgr*2, 1+cgr*2, 1]


def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    d1.targetVol = 3.5 + random.uniform(0.0,0.5)
    d2.targetVol = 3.5 + random.uniform(0.0,0.5)

