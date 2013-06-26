import random
from CellModeller.Signalling.GridDiffusion import GridDiffusion
from CellModeller.Integration.ScipyODEIntegrator import ScipyODEIntegrator
from CellModeller.Integration.CrankNicIntegrator import CrankNicIntegrator
from CellModeller.Integration.CLCrankNicIntegrator import CLCrankNicIntegrator
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy

max_cells = 2**13

def setup(sim):
    dim = (32,32,4)
    size = (5,5,5)
    orig = (-0.5*(dim[0]-1)*size[0],-0.5*(dim[1]-1)*size[1],-0.5*(dim[2]-1)*size[2])
    # Set biophysics, signalling, and regulation models
    sig = GridDiffusion(sim, 2, dim, size, orig, [1e-1,1e-1], [1.0,1.0])
    #integ = ScipyODEIntegrator(sim, 0, 3, max_cells, sig, True)
    integ = CLCrankNicIntegrator(sim, 2, 3, max_cells, sig)
    #integ = CrankNicIntegrator(sim, 2, 3, max_cells, sig)

    biophys = CLBacterium(sim, max_cells=max_cells, max_contacts=32, max_sqs=64**2)
    biophys.addPlane((0,0,-0.5), (0,0,1), 0.75)
    biophys.addPlane((0,0,0.5), (0,0,-1), 1e-2)

    regul = ModuleRegulator(sim, __file__)	# use this file for reg too

    sim.init(biophys, regul, sig, integ)

    sim.addCell(cellType=0, pos=(-20,1,0))
    sim.addCell(cellType=1, pos=(20,0,1))

    # Add some objects to draw the models
    sigrenderer = Renderers.GLGridRenderer(sig, integ)
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(sigrenderer)
    sim.addRenderer(therenderer)

def init(cell):
    #cell.signals[:] = 0.0
    #cell.species[:] = 0.0
    cell.signals = [0,0]
    cell.targetVol = 2.5 + random.uniform(0.0,0.5)
    cell.growthRate = 0.5

def numSignals():
    return 2

def numSpecies():
    return 3

def signalRates(cell, speciesLevels, signalLevels):
    return [0.1, 0.1]

def speciesRates(cell, speciesLevels, signalLevels):
    return [1,1,1]

def update(cells):
    for (id, cell) in cells.items():
        if cell.cellType==0:
            if cell.signals[0]>0.3:
                cell.color = [0.0, 0.0, 1.0]
            else:
                cell.color = [0.1,0.7, 0.1]
        if cell.cellType==1:
            if cell.signals[1]>0.3:
                cell.color = [1.0, 1.0, 0.0]
            else:
                cell.color = [0.7,0.1, 0.1]
        #max(cell.startVol*2.0,0.0): #cell.startvol*2: #random.uniform(1.75,2.0): 
        if cell.volume > cell.targetVol:
            a = 1#random.uniform(0.95,1.05)
            cell.asymm = [a,1]
            cell.divideFlag = True

def divide(parent, d1, d2):
    #for s in range(len(parent.signals)):
    #    u1 = 1 #d1.volume/(d1.volume+d2.volume) #0.5 #+ random.uniform(0,0.01)
    #    u2 = 1 # 1-u1
    #    d1.signals[s] = parent.signals[s] * u1
    #    d2.signals[s] = parent.signals[s] * u2
    d1.targetVol = 2.5 + random.uniform(0.0,0.5)
    d2.targetVol = 2.5 + random.uniform(0.0,0.5)

