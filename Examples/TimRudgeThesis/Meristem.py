from CellModeller.Signalling.GridDiffusion import GridDiffusion
from CellModeller.Integration.CLCrankNicIntegrator import CLCrankNicIntegrator
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import random


grid_dim = (64, 64, 12)
grid_size = (4, 4, 4)
grid_orig = (-128, -128, -8)


def setup(sim):
    biophys = CLBacterium(sim, jitter_z=False)
    reg = ModuleRegulator(sim, sim.moduleName)

    sim.init(biophys, reg, None, None)

    therend = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therend)

    sim.addCell(cellType=0, pos=(0,0,0), len=2.0)


def init(cell):
    cell.target_volume = 3.0
    cell.length = 2.0
    cell.signals = [0]
    cell.species = [0, 0, 0, 0, 0]
    cell.growthRate = 0.0

def max_x_coord(cells):
    mx = -1000
    for i,cell in list(cells.items()):
        mx = max(mx, cell.pos[0])
    return mx

growth_rad = 10.0

def update(cells):
    global growth_rad
    rightmost = -1
    edge = max_x_coord(cells)
    for (i,cell) in list(cells.items()):
        dist = edge - cell.pos[0]
        if dist < growth_rad:
            cell.growthRate = 0.5 #(growth_rad-dist)/growth_rad * 0.5
        else:
            cell.growthRate = 0.0
        if cell.volume > getattr(cell, 'target_volume', 3.0):
            cell.asymm = [1,1]
            cell.divideFlag = True
        cell.color = [cell.growthRate*2*0.9+0.1]*2+[0.5]


def divide(parent, d1, d2):
    d1.target_volume = 3.0 + random.uniform(-0.5, 0.5)
    d2.target_volume = 3.0 + random.uniform(-0.5, 0.5)

def numSignals():
    return 0

def numSpecies():
    return 0
