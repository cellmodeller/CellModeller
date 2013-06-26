from CellModeller.Signalling.GridDiffusion import GridDiffusion
from CellModeller.Integration.CLCrankNicIntegrator import CLCrankNicIntegrator
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
import numpy
import random
import math

max_cells = 2**15

grid_dim = (64, 64, 12)
grid_size = (4, 4, 4)
grid_orig = (-128, -128, -8)


def setup(sim):
    sig = GridDiffusion(sim, 1, grid_dim, grid_size, grid_orig, [1.0])
    integ = CLCrankNicIntegrator(sim, 1, 5, max_cells, sig)
    biophys = CLBacterium(sim, max_cells=max_cells, max_contacts=32, max_sqs=64**2, jitter_z=False, reg_param=2)
    reg = ModuleRegulator(sim, __file__)

    sim.init(biophys, reg, sig, integ)

    sigrend = Renderers.GLGridRenderer(sig, integ)
    therend = Renderers.GL2DBacteriumRenderer(sim)
    sim.addRenderer(sigrend)
    sim.addRenderer(therend)

    sim.addCell(pos=(0,0,0), len=2.0)


def init(cell):
    cell.target_volume = 3.0
    cell.length = 2.0
    cell.signals[:] = [0]
    cell.species[:] = [1000, 0, 0, 0, 0]
    cell.growthRate = 0.5

def numSignals():
    return 1

def numSpecies():
    return 5

# parameters from Danino et al, 2010 SI
cl_preamble = """
const float CA = 1.0;
const float CI = 4.0;//4.0;
const float delta = 1e-3;//delta = 1e-3;
const float alpha = 2500.0;//const float alpha = 2500.0;
const float k = 4.0;//1.0;
const float k1 = 0.1;
const float b = 0.6;
const float gA = 15.f;//15.0;
const float gI = 30.f;//24.0;
const float gH = 0.1f;//0.1;
const float f = 0.3;
const float g = 0.01;
const float D = 2.5f;//2.5;
const float mu = 1.0;
const float mat = 0.005f;
const float gPm = 0.5f;
float A = species[0];
float I = species[1];
float Hi = species[2];
float P = species[3];
float Pm = species[4];
float He = signals[0];
"""

cl_spec_rate = """
rates[0] = CA*Pm - gA*A/(1.f+f*(A+I));
rates[1] = CI*Pm - gI*I/(1.f+f*(A+I));
//rates[2] = b*I/(1.f+k*I) - gH*A*Hi/(1.f+g*A) + D*(He-Hi)*area/volume;
rates[2] = b*I/(1.f+k*I) - gH*A*Hi/(1.f+g*Hi) + D*(He-Hi)*area/volume;
rates[3] = delta + (alpha*Hi*Hi)/(1.f+k1*Hi*Hi) - mat*P;
rates[4] = mat*P - gPm*Pm;
"""
print cl_spec_rate

cl_sig_rate = """
rates[0] = 0;//D*(Hi-He)*area/gridVolume;
"""


def specRateCL():
    return cl_preamble+cl_spec_rate


def sigRateCL():
    return cl_preamble+cl_sig_rate


def update(cells):
    if len(cells) > max_cells-64:
        print 'reached cell limit'
        exit()
    for (i,cell) in cells.items():
        cell.color = [(math.log10(cell.species[0]+1e-12)+6)*0.1, (math.log10(cell.species[1]+1e-12)+6)*0.1, math.log10(cell.species[2]+1e-12)+3]
        if cell.volume > getattr(cell, 'target_volume', 3.0):
            cell.asymm = [1,1]
            cell.divideFlag = True


def divide(parent, d1, d2):
    d1.target_volume = 3.0 + random.uniform(-0.5, 0.5)
    d2.target_volume = 3.0 + random.uniform(-0.5, 0.5)
