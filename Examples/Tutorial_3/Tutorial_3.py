import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import math

from CellModeller.Signalling.GridDiffusion import GridDiffusion #add
from CellModeller.Integration.CLCrankNicIntegrator import CLCrankNicIntegrator #add


max_cells = 10000

#Specify parameter for solving diffusion dynamics #Add
grid_dim = (80, 80, 8) # dimension of diffusion space, unit = number of grid
grid_size = (4, 4, 4) # grid size
grid_orig = (-160, -160, -16) # where to place the diffusion space onto simulation space

n_signals = 2
n_species = 2


def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, jitter_z=False)
    sig = GridDiffusion(sim, n_signals, grid_dim, grid_size, grid_orig, [10.0, 10.0])
    integ = CLCrankNicIntegrator(sim, n_signals, n_species, max_cells, sig)
    
    # use this file for reg too
    regul = ModuleRegulator(sim, sim.moduleName)
    # Only biophys and regulation
    sim.init(biophys, regul, sig, integ)
    
    # Specify the initial cell and its location in the simulation
    sim.addCell(cellType=0, pos=(-3.0,0,0))  #Add
    sim.addCell(cellType=1, pos=(3.0,0,0)) #Add
    
    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)
    sigrend = Renderers.GLGridRenderer(sig, integ) # Add
    sim.addRenderer(sigrend) #Add
    
    sim.pickleSteps = 2

def init(cell):
    # Specify mean and distribution of initial cell size
    cell.targetVol = 3.5 + random.uniform(0.0,0.5)
    # Specify growth rate of cells
    cell.growthRate = 1.0
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
# Db = diffusion rate of beta through the cell membrane

def specRateCL(): # Add if/else, new species
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


def update(cells):
    #Iterate through each cell and flag cells that reach target size for division
    v_max = 0.9
    Km = 0.1 #This can be though of as the degree of mutualism - turns mutualism off when set to 0.
    for (id, cell) in cells.items():
        cell.color = [1,1,1] #[0.1+cell.species[0]/3.0, 0.1+cell.species[1]/3.0, 0.1]
        if cell.cellType==0:
            cell.growthRate = 0.1 + v_max * cell.species[1] / (Km + cell.species[1])
        else:
            cell.growthRate = 0.1 + v_max * cell.species[0] / (Km + cell.species[0])
        
        if cell.volume > cell.targetVol:
            cell.divideFlag = True

def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    d1.targetVol = 3.5 + random.uniform(0.0,0.5)
    d2.targetVol = 3.5 + random.uniform(0.0,0.5)

