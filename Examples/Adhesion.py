import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import math



max_cells = 800
cell_col = [[1.0,0.0,0.0],[0.0,1.0,0.0]]

radius = 0.75/2
delta = 1.9
delta_sig=0.45

adh=100000

def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, max_substeps=10, jitter_z=False, cgs_tol=1e-8, max_contacts=32, max_sqs=192**2,reg_param=2.0, gamma=10.0, max_cells=max_cells, adh_strength=adh)
    
    # use this file for reg too
    regul = ModuleRegulator(sim, sim.moduleName)
    # Only biophys and regulation
    sim.init(biophys, regul, None, None)
    
    # Specify the initial cell and its location in the simulation
    sim.addCell(cellType=0, pos=(0.0,0,0),cellAdh=1.0,rad=radius, length = delta) #Add
    
    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)
    #sigrend = Renderers.GLGridRenderer(sig, integ) # Add
    #sim.addRenderer(sigrend) #Add
    
    sim.pickleSteps = 100

#adhesion logic - here the adhesion strength is the weakest cellAdh of the two contacting cells
def adhLogicCL():
    return """return min(adh_str1, adh_str2);"""

def init(cell):
    # Specify mean and distribution of initial cell size
    #cell.targetVol = 3.5 + random.uniform(0.0,0.5)
    #MODEL 2 constant adder
    cell.targetVol = delta + random.gauss(delta,delta_sig)
    # Specify growth rate of cells
    cell.growthRate = 1.0
    cell.color = cell_col[cell.cellType]

def update(cells):
    #Iterate through each cell and flag cells that reach target size for division
    for (id, cell) in cells.iteritems():
        if cell.volume > cell.targetVol:
            cell.divideFlag = True

def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    #MODEL 1
    #d1.targetVol = rodlen + random.uniform(0.0,1.0)
    #d2.targetVol = rodlen + random.uniform(0.0,1.0)
    #MODEL 2 constant adder
    d1.targetVol = d1.length + random.gauss(delta,delta_sig)
    d2.targetVol = d2.length + random.gauss(delta,delta_sig)

