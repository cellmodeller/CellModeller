from CellModeller.Signalling.GridDiffusion import GridDiffusion
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.Integration.CLCrankNicIntegrator import CLCrankNicIntegrator
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.GUI import Renderers
import numpy
import random



max_cells = 2**11

grid_dim = (64, 64, 12)
grid_size = (4, 4, 4)
grid_orig = (-128, -128, -8)

n_signals = 2
n_species = 5


def setup(sim):
    sig = GridDiffusion(sim, n_signals, grid_dim, grid_size, grid_orig, [0.1, 1.0])
    integ = CLCrankNicIntegrator(sim, n_signals, n_species, max_cells, sig, greensThreshold=1e-13)
    biophys = CLBacterium(sim, max_cells=max_cells, max_contacts=32, max_sqs=32**2, jitter_z=False, reg_param=2)
    reg = ModuleRegulator(sim, __file__)

    sim.init(biophys, reg, sig, integ)

    sigrend = Renderers.GLGridRenderer(sig, integ)
    therend = Renderers.GL2DBacteriumRenderer(sim)
    sim.addRenderer(sigrend)
    sim.addRenderer(therend)

    sim.addCell(pos=(0.0,0,0), len=2.0)


def init(cell):
    cell.target_volume = 3.0
    cell.length = 2.0
    cell.signals = [0.0]*n_species
    cell.species = [0.0]*n_signals
    cell.growthRate = 0.5


def numSignals():
    global n_signals
    return n_signals


def numSpecies():
    global n_species
    return n_species


cl_prefix = \
    """
    // A synthesis
    const float kas1 = 0.01f;
    const float kas2 = 0.1f;

    // I synthesis
    const float kis1 = 0.01f;
    const float kis2 = 0.1f;

    // A production of As/Is
    const float kaa1 = 0.02f;
    const float kaa2 = 0.1f;
    const float kai1 = 0.01f;
    const float kai2 = 0.1f;

    // I production of D
    const float ki1 = 0.1f;
    const float ki2 = 0.01f;

    // D degradation of A/I
    const float kda1 = 0.05f;
    const float kda2 = 0.1f;
    const float kdi1 = 0.05f;
    const float kdi2 = 0.1f;

    // A/I diffusion across the membrane
    const float DA = 0.5f;
    const float DI = 0.5f;

    // basal expression of As,Is,D
    const float eas = 0.001f;
    const float eis = 0.0001f;
    const float ed = 0.00001f;

    // degradation of As, Is, D
    const float das = 0.1f;
    const float dis = 0.1f;
    const float dd = 0.1f;

    //degradation of A, I
    const float da = 0.01f;
    const float di = 0.1f;

    float Ai = species[0];
    float Ii = species[1];
    float As = species[2];
    float Is = species[3];
    float D = species[4];

    float A = signals[0];
    float I = signals[1];
"""


def specRateCL():
    global cl_prefix
    return cl_prefix + """
    rates[0] = kas1*As/(kas2+As) - kda1*D*Ai/(kda2+Ai) + DA*(A-Ai);
    rates[1] = kis1*Is/(kis2+Is) - kdi1*D*Ii/(kdi2+Ii) + DI*(I-Ii);
    rates[2] = eas + kaa1*A*A/(kaa2+A*A) - das*As;
    rates[3] = eis + kai1*A*A/(kai2+A*A) - dis*Is;
    rates[4] = ed + ki1*I*I/(ki2+I*I) - dd*D;
    """

def sigRateCL():
    global cl_prefix
    return cl_prefix + """
    rates[0] = -DA*(A-Ai)*volume - da*A;
    rates[1] = -DI*(I-Ii)*volume - di*I;
    """




def update(cells):
    for i,state in cells.items():
        Ai,Ii,As,Is,D = state.species
        state.color = [D, Ai, Ii]
        if state.volume > getattr(state, 'target_volume', 3.0):
            state.asymm = [1, 1]
            state.divideFlag = True


def divide(parent, d1, d2):
    d1.target_volume = 3.0 + random.uniform(-0.5, 0.5)
    d2.target_volume = 3.0 + random.uniform(-0.5, 0.5)

