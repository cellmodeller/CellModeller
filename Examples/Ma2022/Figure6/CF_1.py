import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import math

from CellModeller.Signalling.GridDiffusion import GridDiffusion #add
from CellModeller.Integration.CLCrankNicIntegrator import CLCrankNicIntegrator #add

delta = 2.0
delta_sig = 0.45
max_cells = 20000
#Specify parameter for solving diffusion dynamics #Add
grid_dim = (120, 120, 12) # dimension of diffusion space, unit = number of grid
grid_size = (4, 4, 4) # grid size
grid_orig = (-240, -240, -24) # where to place the diffusion space onto simulation space

n_signals = 2
n_species = 2 ## Two-way cross-feeding

def setup(sim):
    # Set biophysics module
    biophys = CLBacterium(sim, jitter_z=False, max_cells=max_cells, cgs_tol=1E-5,compNeighbours=True)
    sig = GridDiffusion(sim, n_signals, grid_dim, grid_size, grid_orig, [10.0, 10.0])
    integ = CLCrankNicIntegrator(sim, n_signals, n_species, max_cells, sig)

    # Set up regulation module
    regul = ModuleRegulator(sim, sim.moduleName)
    # Only biophys and regulation
    sim.init(biophys, regul, sig, integ)


    # Specify the initial cell and its location in the simulation
    dist = 5 #distance between cells on the grid
    radius = 60


    for i in range(-19,20):
        for j in range(-19,20):
            R = numpy.sqrt((i*dist)**2 + (j*dist)**2)
            if R <= radius:
                if (i+j)%2 == 0:
                    sim.addCell(cellType=0, pos=(i*dist,j*dist,0), dir=(random.uniform(-1,1),random.uniform(-1,1),0), length=delta,rad=0.4)
                else:
                    sim.addCell(cellType=1, pos=(i*dist,j*dist,0), dir=(random.uniform(-1,1),random.uniform(-1,1),0), length=delta,rad=0.4)

    # Add some objects to draw the models, comment the following four lines when running in terminal.
    #therenderer = Renderers.GLBacteriumRenderer(sim)
    #sim.addRenderer(therenderer)
    #sigrend = Renderers.GLGridRenderer(sig, integ) # Add
    #sim.addRenderer(sigrend) #Add
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
    # Specify initial concentration of chemical species
    cell.species[:] = [0.0]*n_species
    # Specify initial concentration of signaling molecules
    cell.signals[:] = [0.0]*n_signals

cl_prefix = \
    '''
        const float Da = 1.0f;
        const float Db = 1.0f;
        const float ka = 1.f;
        const float kb = 1.f;

        float  alpha_in = species[0];
        float  alpha = signals[0];

        float beta_in = species[1];
        float beta = signals[1];
        '''

# Da = diffusion rate of alpha through the cell membrane

def specRateCL():
    global cl_prefix
    return cl_prefix + '''
        if (cellType==0){
        rates[0] = ka + Da*(alpha-alpha_in)*area/gridVolume;
        rates[1] = Db*(beta-beta_in)*area/gridVolume;

        } else {
        rates[0] = Da*(alpha-alpha_in)*area/gridVolume;
        rates[1] = kb + Db*(beta-beta_in)*area/gridVolume;
        }
        '''

def sigRateCL(): #Add
    global cl_prefix
    return cl_prefix + '''
        rates[0] = -Da*(alpha-alpha_in)*area/gridVolume;
        rates[1] = -Db*(beta-beta_in)*area/gridVolume;

        '''


def update(cells,sim):
    v_max = 0.9
    Km = 0.1 #This can be thought of as the degree of mutualism - turns mutualism off when set to 0.
    #The module regulator now passes the step number into the update function
    #Iterate through each cell and flag cells that reach target size for division
    for (id, cell) in cells.items():

        if cell.volume > cell.targetVol:
            cell.divideFlag = True

        if sim.stepNum < 830:
            if cell.cellType == 0: #if a cell is an acceptor,
                cell.growthRate = 0.1 + v_max * cell.species[1] / (Km + cell.species[1])

            if cell.cellType == 1:
                cell.growthRate = 0.1 + v_max * cell.species[0] / (Km + cell.species[0])


        if sim.stepNum >= 830:
            if cell.cellType == 0: #if a cell is an acceptor,
                cell.growthRate = 0.1 + v_max * cell.species[1] / (Km + cell.species[1])
                for index in cell.neighbours:#loop through all contacts
                    if cells[index].cellType == 1 : #if donors not transconjugants in contact
                        if random.random() < 0.005: # constant probability per unit time
                        #if random.random() < cells[index].effGrowth/50.0: # probability of infection is proportional to donor growth rate
                            sim.HGTevents += 1
                            cell.cellType = 2 #become transconjugant
                            cell.growthRate = 0.1 + v_max * cell.species[1] / (Km + cell.species[1])
                            cell.color = (0.98, 1.0, 0.15) ## yellow

            if cell.cellType == 1:
                cell.growthRate = 0.1 + v_max * cell.species[0] / (Km + cell.species[0])


def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    d1.targetVol = d1.length + random.gauss(delta,delta_sig)
    d2.targetVol = d2.length + random.gauss(delta,delta_sig)
