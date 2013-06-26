import random
from CellModeller.Signalling.GridDiffusion import GridDiffusion 
from CellModeller.Integration.ScipyODEIntegrator import ScipyODEIntegrator
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CGSBacterium import CGSBacterium
from CellModeller.GUI import Renderers
import numpy


def setup(sim):
    dim = (128,128,1)
    size = (20,20,20)
    orig = (-0.5*dim[0]*size[0],-0.5*dim[1]*size[1],-0.5*dim[2]*size[2])
    # Set biophysics, signalling, and regulation models
    sig = GridDiffusion(sim, 1, dim, size, orig, [1e-1])
    integ = ScipyODEIntegrator(sim, 0, 2, 4096, sig, True)
    biophys = CGSBacterium(sim)
    regul = ModuleRegulator(sim, __file__)	# use this file for reg too

    sim.init(biophys, regul, sig, integ)
    sim.addCell()

    # Add some objects to draw the models
    #sigrenderer = Renderers.GLGridRenderer(sig, integ)
    therenderer = Renderers.GLCelBacteriumRenderer(sim)
    #sim.addRenderer(sigrenderer)
    sim.addRenderer(therenderer)

def init(cell):
    cell.signals = [0]
    cell.species = [0, 0]

def numSignals():
    return 1

def numSpecies():
    return 2

def signalRates(cell, speciesLevels, signalLevels):
    return [1.0]

def speciesRates(cell, speciesLevels, signalLevels):
    return [1.0, 2.0]

def update(cells):
    for (id, cell) in cells.items():
        if cell.volume > max(cell.startVol*2.0,0.0): #cell.startvol*2: #random.uniform(1.75,2.0):
            a = 1#random.uniform(0.95,1.05)
            cell.asymm = [a,1]
            cell.divideFlag = True

def divide(parent, d1, d2):
    for s in range(len(parent.signals)):
        u1 = 1 #d1.volume/(d1.volume+d2.volume) #0.5 #+ random.uniform(0,0.01)
        u2 = 1 # 1-u1
        d1.signals[s] = parent.signals[s] * u1
        d2.signals[s] = parent.signals[s] * u2

