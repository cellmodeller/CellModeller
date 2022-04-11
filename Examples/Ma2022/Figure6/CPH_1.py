import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import math

delta = 2.0
delta_sig = 0.45
max_cells = 20000

def setup(sim):
    # Set biophysics module
    biophys = CLBacterium(sim, jitter_z=False, max_cells=max_cells, cgs_tol=1E-5,compNeighbours=True)

    # Set up regulation module
    regul = ModuleRegulator(sim, sim.moduleName)
    # Only biophys and regulation
    sim.init(biophys, regul, None, None)

    # Specify the initial cell and its location in the simulation
    dist = 5 #distance between cells on the grid
    radius = 60


    for i in range(-19,20):
        for j in range(-19,20):
            R = numpy.sqrt((i*dist)**2 + (j*dist)**2)
            if R <= radius:
                if (i+j)%2==0:
                    sim.addCell(cellType=0, pos=(i*dist,j*dist,0), dir=(random.uniform(-1,1),random.uniform(-1,1),0), length=delta,rad=0.4)
                else:
                    sim.addCell(cellType=1, pos=(i*dist,j*dist,0), dir=(random.uniform(-1,1),random.uniform(-1,1),0), length=delta,rad=0.4) ## inoculate recipients randomly.

            #elif R < 18:
            #    sim.addCell(cellType=1, pos=(i*dist,j*dist,0), dir=(random.uniform(-1,1),random.uniform(-1,1),0), length=delta,rad=0.4)


     #producer at frontier: in total 56 donors; and 249 recipients, to make it roughly 1:1, need to add 97 more donors.



    #ACCEPTORS:
    #sim.addCell(cellType=0, pos=(-5,-5,0), dir=(1,0,0), length=delta,rad=0.4)
    #sim.addCell(cellType=0, pos=(5,5,0), dir=(1,0,0), length=delta,rad=0.4)
    #DONORS:
    #sim.addCell(cellType=1, pos=(5,-5,0), dir=(1,0,0), length=delta,rad=0.4)
    #sim.addCell(cellType=1, pos=(-5,5,0), dir=(1,0,0), length=delta,rad=0.4)

    # Add some objects to draw the models
    #therenderer = Renderers.GLBacteriumRenderer(sim)
    #sim.addRenderer(therenderer)
    #sigrend = Renderers.GLGridRenderer(sig, integ) # Add
    #sim.addRenderer(sigrend) #Add

    # Add some objects to draw the models
    #therenderer = Renderers.GLBacteriumRenderer(sim)
    #sim.addRenderer(therenderer)
    # Specify how often data is saved
    sim.pickleSteps = 10

    sim.countHGT = True #Turns on HGT event counting
    sim.HGTevents = 0 #You need to start the counter here



def init(cell):
    # Specify mean and distribution of initial cell size
    cell.targetVol = cell.length + random.gauss(delta,delta_sig)
    # Specify growth rate of cells
    cell.growthRate = 1.0
    #if cell.cellType == 0:s
    #    cell.color = (1.0, 0.1, 0.1) #red acceptor, grows faster
    #    cell.growthRate = 1.1
    if cell.cellType == 0:
        cell.color = (1.0, 0.1, 0.1)#red consumer, acceptor
    elif cell.cellType == 1:
        cell.color = (0.32, 1.0, 0.96) #cyan producer, donor
    elif cell.cellType == 2:
        cell.color = (0.98, 1.0, 0.15) #yellow transconjugant

def update(cells,sim): #The module regulator now passes the step number into the update function
    #Iterate through each cell and flag cells that reach target size for division
    for (id, cell) in cells.items():
        if cell.volume > cell.targetVol:
            cell.divideFlag = True


        if sim.stepNum < 400:
            cell.growthRate = 1.0

        if sim.stepNum >= 400: ### Apply HGT starting at time step 160
            if cell.cellType == 0: #if a cell is an acceptor,
                for index in cell.neighbours:#loop through all contacts
                    if cells[index].cellType == 1 : #if donors not transconjugants in contact
                        if random.random() < 0.005: # constant probability per unit time
                        #if random.random() < cells[index].effGrowth/50.0: # probability of infection is proportional to donor growth rate
                            sim.HGTevents += 1
                            cell.cellType = 2 #become transconjugant
                            cell.growthRate = 1.0
                            cell.color = (0.98, 1.0, 0.15)


def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    d1.targetVol = d1.length + random.gauss(delta,delta_sig)
    d2.targetVol = d2.length + random.gauss(delta,delta_sig)
