import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import math








max_cells = 400000

def setup(sim):
    # Set biophysics, signalling, and regulation models
    # Ton: [KNOB] jitter_z = True would allow cells to pile up on the z-axis, jitter_z = False would keep everything on 2Dplane
    biophys = CLBacterium(sim, max_substeps=8, max_cells=max_cells, max_contacts=32, max_sqs=192**2, jitter_z=False, reg_param=2, gamma=10)
 


    regul = ModuleRegulator(sim, __file__)	# use this file for reg too
    # Only biophys and regulation
    sim.init(biophys, regul, None, None)

    
    #Ton: [KNOB] for changing initial location and type of cell
    sim.addCell(cellType=0, pos=(-8,0,0)) 
    sim.addCell(cellType=1, pos=(8,0,0)) 
    #Ton: when do I specify len?
    #Ton: can I also set up orientation of initial cell?

    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)

    #Ton: [KNOB?] for chaning time se
    sim.pickleSteps = 10  #Ton tried changing this but printed out time step does not change. why?


def init(cell):
    #Ton: [KNOB] initial cell size
    cell.targetVol = 2.5 + random.uniform(0.0,0.5)
    #Ton: [KNOB] cell growth rate
    cell.growthRate = 2.0

def numSignals():
    return 0

def numSpecies():
    return 0

def update(cells):
    for (id, cell) in cells.iteritems():
        cell.color = [cell.cellType*0.6+0.1, 1.0-cell.cellType*0.6, 0.3]
        #max(cell.startVol*2.0,0.0): #cell.startvol*2: #random.uniform(1.75,2.0):
        if cell.volume > cell.targetVol:
            a = 1#random.uniform(0.95,1.05)
            cell.asymm = [a,1]
            cell.divideFlag = True

def divide(parent, d1, d2):
   #Ton: [KNOB] target cell size that triggers cell division
    d1.targetVol = 2.5 + random.uniform(0.0,0.5)
    d2.targetVol = 2.5 + random.uniform(0.0,0.5)

