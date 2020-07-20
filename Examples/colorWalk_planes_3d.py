import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
import numpy
import math

max_cells = 10000
rodlen = 3.5

def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, max_planes=9)

    sim.dt = 0.025
    #biophys.addPlane((0,0,-0.5), (0,0,1), 1.0)
    #biophys.addPlane((0,0,0.5), (0,0,-1), math.sqrt(7.5e-4))

    '''
    angs = [i*math.pi/8 for i in range(8)]
    for a in angs:
        x = math.cos(a)
        y = math.sin(a)
        p = (30*x,30*y,0)
        n = (-x,-y,0)
        biophys.addPlane(p,n, 1.0)
    '''
    biophys.addPlane((0,0,0),(0,0,-1),1.0)

    regul = ModuleRegulator(sim)	# use this file for reg too
    # Only biophys and regulation
    sim.init(biophys, regul, None, None)

    #ct = 0
    #for x in range(-4,5):
    #    for y in range(-4,5):
    sim.addCell(cellType=0, pos=(0,0,-5), color=(1,1,1)) #x*40,y*40,0))
    #        ct += 1
    #sim.addCell(cellType=1, pos=(20,0,0))


    if sim.is_gui:
        # Add some objects to draw the models
        from CellModeller.GUI import Renderers
        therenderer = Renderers.GLBacteriumRenderer(sim)
        sim.addRenderer(therenderer)

    sim.pickleSteps = 10

def init(cell):
    cell.targetVol = rodlen + random.uniform(0.0,0.5)
    cell.growthRate = 1.0

def numSignals():
    return 0

def numSpecies():
    return 0

def update(cells):
    for (id, cell) in cells.items():
        if cell.volume > cell.targetVol:
            cell.divideFlag = True

def divide(parent, d1, d2):
    d1.targetVol = rodlen + random.uniform(0.0,0.5)
    d2.targetVol = rodlen + random.uniform(0.0,0.5)
    u1 = numpy.random.uniform(-1.0,1.0,size=(3,))
    u2 = numpy.ones((3,)) - u1
    d1.color = [parent.color[i] + 0.1*u1[i] for i in range(3)]

