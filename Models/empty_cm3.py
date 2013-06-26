import random 
from CellModeller.Integration.ScipyODEIntegrator import ScipyODEIntegrator
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.PlantModels.CM3 import CM3
from CellModeller.GUI import Renderers
import numpy


def setup(sim):
    # Set biophysics, signalling, and regulation models
    integ = ScipyODEIntegrator(sim, 0, 1, 4096)
    biophys = CM3(sim)
    regul = ModuleRegulator(sim, __file__)	# use this file for reg too

    sim.init(biophys, regul, None, integ)
    sim.addCell()
    
    # Add some objects to draw the models
    therenderer = Renderers.GLPlantRenderer(sim)
    sim.addRenderer(therenderer)

def init(cell):
    cell.signals = []
    cell.species = [0.0]
    cell.targetVol = 400.0
    
def numSignals():
    return 0  

def numSpecies():
    return 1

def signalRates(cell, speciesLevels, signalLevels):
    return []

def speciesRates(cell, speciesLevels, signalLevels):
    return [1.0]
    
def update(cells):
    for (id, cell) in cells.items():
        if cell.volume > cell.targetVol:
            a = 1#random.uniform(0.95,1.05)
            cell.asymm = [a,1]
            cell.divideFlag = True

def divide(p, d1, d2):
    d1.targetVol = 200.0+random.uniform(0.0,200.0)
    d2.targetVol = 200.0+random.uniform(0.0,200.0)
    
