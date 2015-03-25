import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import math

from CellModeller.Signalling.GridDiffusion import GridDiffusion #add
from CellModeller.Integration.CLCrankNicIntegrator import CLCrankNicIntegrator #add
from CellModeller.Integration.ScipyODEIntegrator import ScipyODEIntegrator #add


max_cells = 400000
#Specify parameter for solving diffusion dynamics
grid_dim = (64, 8, 12) # dimension of diffusion space, unit = number of grid
grid_size = (4, 4, 4) # grid size
grid_orig = (-128, -14, -8) # where to place the diffusion space onto simulation space


def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, max_substeps=8, max_cells=max_cells, max_contacts=32, max_sqs=192**2, jitter_z=False, reg_param=2, gamma=10)
 
    # Ton: [KNOB] add the planes to set physical  boundaries of cell growth
    biophys.addPlane((0,-16,0), (0,1,0), 1)
    biophys.addPlane((0,16,0), (0,-1,0), 1)

    sig = GridDiffusion(sim, 1, grid_dim, grid_size, grid_orig, [10.0])
    integ = CLCrankNicIntegrator(sim, 1, 5, max_cells, sig, boundcond='reflect')

    # use this file for reg too
    regul = ModuleRegulator(sim, __file__)	
    # Only biophys and regulation
    sim.init(biophys, regul, sig, integ)


    # Specify the initial cell and its location in the simulation
    sim.addCell(cellType=1, pos=(-20.0,0,0), len=2.0) #Add
    sim.addCell(cellType=0, pos=(20.0,0,0), len=2.0) #Add



    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)
    sigrend = Renderers.GLGridRenderer(sig, integ) # Add
    sim.addRenderer(sigrend) #Add

    sim.pickleSteps = 10



def init(cell):
    # Specify mean and distribution of initial cell size
    cell.targetVol = 2.5 + random.uniform(0.0,0.5)
    # Specify growth rate of cells
    cell.growthRate = 2.0

def numSignals(): # We leave this part empty as we do not simulate cell-cell signaling here
    return 1

def numSpecies(): # We leave this part empty as we do not simulate intracellular gene network here
    return 5

def specRateCL():
    return '''
    const float D1 = 0.1f;
    const float d1 = 1e-3;
    const float k1 = 1.f;
    const float k2 = 1.f;
    const float g1 = 1e-1;
    const float d2 = 1e-5;
    const float k3 = 1.f;
    const float k4 = 5e-5;
    const float k5 = 0.1f;
    const float k6 = 1e-2;
    const float d3 = 1e-3;
    const float g2 = 1e-1;
    const float dr = 0.1f;

    float AHL = signals[0];
    float LuxI = species[0];
    float AHLi = species[1];
    float AiiA = species[2];
    float CFP = species[3];
    float YFP = species[4];

    if (cellType==0) {
        rates[0] = 0.f;
        rates[1] = k1*LuxI/(k2+LuxI) - k5*AiiA*AHLi/(k6+AHLi) + D1*(AHL-AHLi)*area/volume;
        rates[2] = d3 - g2*AiiA;
        rates[3] = d2 + k3*AHLi*AHLi/(k4+AHLi*AHLi) - dr*CFP;
        rates[4] = 0.f;
    } else {
        rates[0] = d1 - g1*LuxI;
        rates[1] = k1*LuxI/(k2+LuxI) + D1*(AHL-AHLi)*area/volume;
        rates[2] = 0.f;
        rates[3] = 0.f;
        rates[4] = d2 + k3*AHLi*AHLi/(k4+AHLi*AHLi) - dr*YFP;
    }
    '''

def sigRateCL():
    return '''
    const float D1=0.1f;
    float AHL = signals[0];
    float AHLi = species[1];
    rates[0] = -D1*(AHL-AHLi)*area/gridVolume;
    '''





def update(cells):
    #Iterate through each cell and flag cells that reach target size for division
    for (id, cell) in cells.iteritems():
        cell.color = [cell.cellType*0.6+0.1, 1.0-cell.cellType*0.6, 0.3]
        if cell.volume > cell.targetVol:
            a = 1
            cell.asymm = [a,1]
            cell.divideFlag = True

def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    d1.targetVol = 2.5 + random.uniform(0.0,0.5)
    d2.targetVol = 2.5 + random.uniform(0.0,0.5)

