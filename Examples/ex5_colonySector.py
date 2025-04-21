import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.Biophysics.BacterialModels.Bacterium import Bacterium
from CellModeller.GUI import Renderers
import numpy
import math

N0 = 10
colors = {0: [1, 0, 0], 1: [0, 1, 0]}
radius = {0: 0.5, 1: 1}

def setup(sim):
    sim.dt = 0.025

    # Set biophysics, signalling, and regulation models
    #biophys = CLBacterium(sim, jitter_z=False, gamma = 1, max_cells=100000)
    biophys = Bacterium(sim, gamma_factor=5, muA=2, sub_steps=3)

    regul = ModuleRegulator(sim, sim.moduleName)	# use this file for reg too
    # Only biophys and regulation
    sim.init(biophys, regul, None, None)

    #biophys.addPlane((0,0,0),(0,0,1),1.0) #Base plane  
    #biophys.addPlane((10,0,0),(-1,0,0),1.0)
    #biophys.addPlane((-10,0,0),(1,0,0),1.0)
    #biophys.addPlane((0,10,0),(0,-1,0),1.0)
    #biophys.addPlane((0,-10,0),(0,1,0),1.0)

    for _ in range(1):
        ct = 0 # (random.uniform(0, 1) > 0.5) * 1
        R = 0 # random.uniform(0, 100)
        theta = random.uniform(0, 2 * numpy.pi)
        pos = R * numpy.array([numpy.cos(theta), numpy.sin(theta), 0])
        dir = numpy.random.uniform(0, 100, size=(3,))
        dir[2] = 0
        dir = dir / numpy.linalg.norm(dir)
        sim.addCell(cellType=ct, pos=tuple(pos), dir=tuple(dir), rad=radius[ct])
    
    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)
    sim.pickleSteps = 1

def init(cell):
    cell.targetVol = 5.5 + random.uniform(0.0,0.5) if cell.cellType==1 else 3.5 + random.uniform(0.0,0.5)
    cell.growthRate = 1.0
    cell.n_a = N0//2
    cell.n_b = N0 - cell.n_a

def update(cells):
    for (id, cell) in cells.items():
        #cell.color = [0.1, cell.n_a/3.0, cell.n_b/3.0]
        cell.color = colors[cell.cellType]
        if cell.volume > cell.targetVol:
            cell.divideFlag = True

def divide(parent, d1, d2):
    d1.targetVol = 5.5 + random.uniform(0.0,0.5) if d1.cellType==1 else 3.5 + random.uniform(0.0,0.5)
    d2.targetVol = 5.5 + random.uniform(0.0,0.5) if d2.cellType==1 else 3.5 + random.uniform(0.0,0.5)
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
