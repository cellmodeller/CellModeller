import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import math

cell_cols = {0:[0,1.0,0], 1:[1.0,0,0], 2:[0,0,1.0]} #RGB cell colours
cell_lens = {0:1.0, 1:2.0, 2:3.5} #target cell lengths
cell_growr = {0:2.0, 1:1.1, 2:0.8} #growth rates

def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, jitter_z=False)

    # use this file for reg too
    regul = ModuleRegulator(sim, sim.moduleName)
    # Only biophys and regulation
    sim.init(biophys, regul, None, None)
 
    # Specify the initial cell and its location in the simulation
    sim.addCell(cellType=0, pos=(0,0,0), dir=(1,0,0))
    sim.addCell(cellType=1, pos=(6,0,0), dir=(1,0,0))
    sim.addCell(cellType=2, pos=(-6,0,0), dir=(1,0,0))


    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)
    sim.pickleSteps = 20

def max_y_coord(cells):
    #finds the largest y-coordinate in the colony
    my = 0.0
    for i,cell in list(cells.items()):
        my = max(my, cell.pos[1])
    return my

def init(cell):
    # Specify mean and distribution of initial cell size
    cell.targetVol = cell_lens[cell.cellType] + random.uniform(0.0,0.5)
    # Specify growth rate of cells
    cell.growthRate = cell_growr[cell.cellType]
    cell.color = cell_cols[cell.cellType]

def update(cells):
    #Iterate through each cell and flag cells that reach target size for division
    maxy = max_y_coord(cells)
    for (id, cell) in cells.items():
        dist = maxy - cell.pos[1]
        growthZone = 5.0 #width of the growth zone
        if dist < growthZone:
            cell.growthRate = cell_growr[cell.cellType]
        else:
            cell.color = numpy.divide(cell_cols[cell.cellType],2)
            cell.growthRate = 0.0
        if cell.volume > cell.targetVol:
            cell.divideFlag = True

def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    d1.targetVol = cell_lens[parent.cellType] + random.uniform(0.0,0.5)
    d2.targetVol = cell_lens[parent.cellType] + random.uniform(0.0,0.5)