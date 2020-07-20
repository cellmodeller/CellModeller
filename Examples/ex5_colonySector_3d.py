import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
import numpy
import math


def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, jitter_z=True)

    regul = ModuleRegulator(sim, sim.moduleName)	# use this file for reg too
    # Only biophys and regulation
    sim.init(biophys, regul, None, None)

    sim.addCell(cellType=0, pos=(0,0,0))
    biophys.addPlane((0,0,0),(0,0,-1),1.0)
    biophys.addPlane((0,0,0),(0,0,1),5e-2)

    if sim.is_gui:
        # Add some objects to draw the models
        from CellModeller.GUI import Renderers
        therenderer = Renderers.GLBacteriumRenderer(sim)
        sim.addRenderer(therenderer)

    sim.pickleSteps = 10

def init(cell):
    cell.targetVol = 3.5 + random.uniform(0.0,0.5)
    cell.growthRate = 1.0
    cell.n_a = 3
    cell.n_b = 3

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
    for p in plasmids[:6]:
        if p == 0: d1.n_a +=1
        else: d1.n_b +=1
    for p in plasmids[6:12]:
        if p == 0: d2.n_a +=1
        else: d2.n_b +=1
    assert parent.n_a + parent.n_b == 6
    assert d1.n_a + d1.n_b == 6
    assert d2.n_a + d2.n_b == 6
    assert parent.n_a*2 == d1.n_a+d2.n_a
    assert parent.n_b*2 == d1.n_b+d2.n_b
    assert parent.n_a > 0 or (d1.n_a == 0 and d2.n_a == 0)
    assert parent.n_b > 0 or (d1.n_b == 0 and d2.n_b == 0)
