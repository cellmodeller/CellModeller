from CellModeller.Signalling.GridDiffusion import GridDiffusion
from CellModeller.Integration.CLCrankNicIntegrator import CLCrankNicIntegrator
from CellModeller.Integration.ScipyODEIntegrator import ScipyODEIntegrator
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import random

max_cells = 2**15
grid_dim = (64, 64, 12)
grid_size = (4, 4, 4)
grid_orig = (-128, -128, -8)


def setup(sim):
    sig = GridDiffusion(sim, 1, grid_dim, grid_size, grid_orig, [10.0], initLevels=[5e-3])
    #integ = ScipyODEIntegrator(sim, 1, 4, max_cells, sig, True)
    integ = CLCrankNicIntegrator(sim, 1, 5, max_cells, sig)
    biophys = CLBacterium(sim, max_cells=max_cells, max_contacts=32, max_sqs=32**2, jitter_z=True, reg_param=2)
    biophys.addPlane((0,0,-0.5), (0,0,1), 0.2)
    biophys.addPlane((0,0,0.5), (0,0,-1), 1e-3)

    reg = ModuleRegulator(sim, __file__)

    sim.init(biophys, reg, sig, integ)

    sigrend = Renderers.GLGridRenderer(sig, integ)
    therend = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(sigrend)
    sim.addRenderer(therend)


    sim.addCell(cellType=0, pos=(0,0,0))

def init(cell):
    cell.target_volume = 4.0
    cell.signals[:] = [0]
    cell.species[:] = [0, 0, 0, 0, 0]
    cell.growthRate = 0.5

def numSignals():
    return 1

def numSpecies():
    return 5



def specRateCL():
    return '''
    const float D1 = 1.f;
    const float d1 = 1e-3;
    const float k1 = 1.f;
    const float k2 = 1.f;
    const float g1 = 1e-1;
    const float d2 = 1e-5;
    const float k3 = 1.f;
    const float k4 = 1e-6; // switching point of Lux promoter**2 
    const float k5 = 2.f;
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
    rates[0] = 0.f;
    rates[1] = - k5*AiiA*AHLi/(k6+AHLi) + D1*(AHL-AHLi)*area/volume;
    rates[2] = d3 - g2*AiiA;
    rates[3] = d2 + k3*AHLi*AHLi/(k4+AHLi*AHLi) - g2*CFP;
    rates[4] = 0.f;
    '''

def sigRateCL():
    return '''
    const float D1=1.f;
    float AHL = signals[0];
    float AHLi = species[1];
    rates[0] = -D1*(AHL-AHLi)*area/gridVolume;
    '''


def update(cells):
    for (i,cell) in cells.items():
        cell.color = [0.1, 0.1+cell.species[3]/10.0, 0.1+cell.species[4]/1e-2]
        if cell.volume > getattr(cell, 'target_volume', 4.0):
            cell.asymm = [1,1]
            cell.divideFlag = True



def divide(parent, d1, d2):
    d1.target_volume = 3.0 + random.uniform(-0.5, 0.5)
    d2.target_volume = 3.0 + random.uniform(-0.5, 0.5)

