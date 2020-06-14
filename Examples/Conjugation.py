import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import math

delta = 1.9
delta_sig = 0.45

def setup(sim):
    # Set biophysics module
    biophys = CLBacterium(sim, jitter_z=False, max_cells=50000,gamma=10.0, cgs_tol=1E-5,compNeighbours=True)

    # Set up regulation module
    regul = ModuleRegulator(sim, sim.moduleName)	
    # Only biophys and regulation
    sim.init(biophys, regul, None, None)
 
    # Specify the initial cell and its location in the simulation
    sim.addCell(cellType=0, pos=(-5,0,0), dir=(1,0,0), length=delta,rad=0.4) #acceptor
    sim.addCell(cellType=1, pos=(5,0,0), dir=(1,0,0), length=delta,rad=0.4) #donor

    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)
    
    # Specify how often data is saved
    sim.pickleSteps = 50

def init(cell):
    # Specify mean and distribution of initial cell size
    cell.targetVol = cell.length + random.gauss(delta,delta_sig)
    # Specify growth rate of cells
    cell.growthRate = 1.0
    if cell.cellType == 0: cell.color = (1.0,0.0,0.0)
    elif cell.cellType == 1: cell.color = (0.0,1.0,0.0)
    elif cell.cellType == 2: cell.color = (0.0,0.0,1.0)


def update(cells):
    #Iterate through each cell and flag cells that reach target size for division
    for (id, cell) in cells.items():
        if cell.volume > cell.targetVol:
            cell.divideFlag = True
        
        infect_chances = 0
        if cell.cellType == 0: #if a cell is an acceptor,
            for index in cell.neighbours:#loop through all contacts
                if cells[index].cellType != 0: #if donor or transconjugant is in contact
                    #if random.random() < 0.01: # constant probability per unit time
                    if random.random() < cells[index].effGrowth/10.0: # probability of infection is proportional to donor growth rate
                        cell.cellType = 2 #become transconjugant

        if cell.cellType == 0: cell.color = (1.0,0.0,0.0)
        elif cell.cellType == 1: cell.color = (0.0,1.0,0.0)
        elif cell.cellType == 2: cell.color = (0.0,0.0,1.0)



def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    d1.targetVol = d1.length + random.gauss(delta,delta_sig)
    d2.targetVol = d2.length + random.gauss(delta,delta_sig)

