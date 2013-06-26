import random 
from CellModeller.Integration.ScipyODEIntegrator import ScipyODEIntegrator
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.PlantModels.CM3 import CM3
from CellModeller.Signalling.NeighbourDiffusion import NeighbourDiffusion
from CellModeller.GUI import Renderers
import numpy


def setup(sim):
    # Set biophysics, signalling, and regulation models
    sig = NeighbourDiffusion(sim, 1, [0.001])
    integ = ScipyODEIntegrator(sim, 1, 1, 4096, sig)
    biophys = CM3(sim)
    regul = ModuleRegulator(sim, __file__)	# use this file for reg too

    sim.init(biophys, regul, sig, integ)
    sim.addCell()
    
    # Add some objects to draw the models
    therenderer = Renderers.GLPlantRenderer(sim)
    sigrenderer = Renderers.GLPlantSignalRenderer(sim, [0])
    sim.addRenderer(sigrenderer)
    sim.addRenderer(therenderer)

def init(cell):
    cell.signals = [0.0]
    cell.species = [0.0]
    
def numSignals():
    return 1  

def numSpecies():
    return 1

def signalRates(cell, speciesLevels, signalLevels):
    if cell.pos[0]>0:
        return [1.0]
    else:
        return [0.0]

def speciesRates(cell, speciesLevels, signalLevels):
    return [1.0]
    
def update(cells):
    for (id, cell) in cells.items():
        if cell.volume > max(cell.startVol*2.0,0.0): #cell.startvol*2: #random.uniform(1.75,2.0):
            a = 1#random.uniform(0.95,1.05)
            cell.asymm = [a,1]
            cell.divideFlag = True
    
