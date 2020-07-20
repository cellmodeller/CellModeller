import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import math

N0 = 10

def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, jitter_z=False, gamma = 100, max_cells=100000, max_planes=1)

    regul = ModuleRegulator(sim, sim.moduleName)	# use this file for reg too
    # Only biophys and regulation
    sim.init(biophys, regul, None, None)

    #biophys.addPlane((0,0,0),(0,0,1),1.0) #Base plane
    #biophys.addPlane((10,0,0),(-1,0,0),1.0)
    #biophys.addPlane((-10,0,0),(1,0,0),1.0)
    #biophys.addPlane((0,10,0),(0,-1,0),1.0)
    #biophys.addPlane((0,-10,0),(0,1,0),1.0)

    sim.addCell(cellType=0, pos=(0,0,0))

    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)
    sim.pickleSteps = 1

def init(cell):
    cell.targetVol = 3.5 + random.uniform(0.0,0.5)
    cell.growthRate = 1.0
    cell.n_a = N0//2
    cell.n_b = N0 - cell.n_a

def update(cells):
    for (id, cell) in cells.items():
        cell.color = [0.1, cell.n_a/3.0, cell.n_b/3.0]
        if cell.volume > cell.targetVol:
            cell.divideFlag = True

def divide(parent, d1, d2):
    d1.targetVol = 3.5 + random.uniform(0.0,0.5)
    d2.targetVol = 3.5 + random.uniform(0.0,0.5)
    plasmids = [0]*parent.n_a*2 + [1]*parent.n_b*2
    random.shuffle(plasmids)
    d1.n_a = 0
    d1.n_b = 0
    d2.n_a = 0
    d2.n_b = 0
    for p in plasmids[:N0]:
        if p == 0: d1.n_a +=1
        else: d1.n_b +=1
    for p in plasmids[N0:2*N0]:
        if p == 0: d2.n_a +=1
        else: d2.n_b +=1
    assert parent.n_a + parent.n_b == N0
    assert d1.n_a + d1.n_b == N0
    assert d2.n_a + d2.n_b == N0
    assert parent.n_a*2 == d1.n_a+d2.n_a
    assert parent.n_b*2 == d1.n_b+d2.n_b
    assert parent.n_a > 0 or (d1.n_a == 0 and d2.n_a == 0)
    assert parent.n_b > 0 or (d1.n_b == 0 and d2.n_b == 0)
