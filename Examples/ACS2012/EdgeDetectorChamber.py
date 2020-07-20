
import random
from CellModeller.Signalling.GridDiffusion import GridDiffusion
from CellModeller.Integration.CLCrankNicIntegrator import CLCrankNicIntegrator
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy


max_cells = 2**15

grid_dim = (64, 8, 12)
grid_size = (4, 4, 4)
grid_orig = (-128, -14, -8)


def setup(sim):
    sig = GridDiffusion(sim, 1, grid_dim, grid_size, grid_orig, [10.0])
    integ = CLCrankNicIntegrator(sim, 1, 5, max_cells, sig, boundcond='reflect')
    biophys = CLBacterium(sim, max_cells=max_cells, jitter_z=False)
    biophys.addPlane((0,-16,0), (0,1,0), 1)
    biophys.addPlane((0,16,0), (0,-1,0), 1)
    reg = ModuleRegulator(sim)

    sim.init(biophys, reg, sig, integ)

    sim.addCell(cellType=1, pos=(-20.0,0,0), len=2.0)
    sim.addCell(cellType=0, pos=(20.0,0,0), len=2.0)

    sigrend = Renderers.GLGridRenderer(sig, integ)
    therend = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(sigrend)
    sim.addRenderer(therend)

    sim.pickleSteps = 10


def init(cell):
    cell.targetVol = 3.5 + random.uniform(0.0,0.5)
    cell.length = 2.0
    cell.signals[:] = [0]
    cell.species[:] = [0, 0, 0, 0, 0]
    cell.growthRate = 1 

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
    if len(cells) > max_cells-16:
        print('reached cell limit')
        exit()
    for (i,cell) in cells.items():
        cell.color = [0.1, 0.1+cell.species[3]/10.0, 0.1+cell.species[4]*20.0]
        if cell.volume > getattr(cell, 'target_volume', 3.0):
            cell.divideFlag = True

def divide(parent, d1, d2):
    d1.targetVol = 3.5 + random.uniform(0.0,0.5)
    d2.targetVol = 3.5 + random.uniform(0.0,0.5)
