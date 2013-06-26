import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import math

max_cells = 400000

def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, max_substeps=8, max_cells=max_cells, max_contacts=32, max_sqs=128**2, jitter_z=False, reg_param=0.04, gamma=25)

    regul = ModuleRegulator(sim, __file__)	# use this file for reg too
    # Only biophys and regulation
    sim.init(biophys, regul, None, None)

    sim.addCell(cellType=0, pos=(0,0,0))

    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)
    sim.pickleSteps = 10

def init(cell):
    cell.targetVol = 3.5 + random.uniform(0.0,0.5)
    cell.growthRate = 1.0
    cell.n_a = 2
    cell.n_b = 2

def numSignals():
    return 0

def numSpecies():
    return 0

def update(cells):
    for (id, cell) in cells.iteritems():
        if cell.n_a > 0 and cell.n_b == 0:
            cell.color = [0.2,0.8,0.2]
        elif cell.n_a == 0 and cell.n_b > 0:
            cell.color = [0.2,0.2,0.8]
        else:
            cell.color = [0.2,0.2,0.2]
        #max(cell.startVol*2.0,0.0): #cell.startvol*2: #random.uniform(1.75,2.0):
        if cell.volume > cell.targetVol:
            a = 1#random.uniform(0.95,1.05)
            cell.asymm = [a,1]
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
    for p in plasmids[:4]:
        if p == 0: d1.n_a +=1
        else: d1.n_b +=1
    for p in plasmids[4:8]:
        if p == 0: d2.n_a +=1
        else: d2.n_b +=1
    assert parent.n_a + parent.n_b == 4
    assert d1.n_a + d1.n_b == 4
    assert d2.n_a + d2.n_b == 4
    assert parent.n_a*2 == d1.n_a+d2.n_a
    assert parent.n_b*2 == d1.n_b+d2.n_b
    assert parent.n_a > 0 or (d1.n_a == 0 and d2.n_a == 0)
    assert parent.n_b > 0 or (d1.n_b == 0 and d2.n_b == 0)

