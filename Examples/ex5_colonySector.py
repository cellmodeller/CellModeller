import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
import numpy as np
import math

l = 60.
K = 10.
n = 1.5
Kn = K**n
_b_ = 10.
_N0_ = 10.

def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, jitter_z=False)

    regul = ModuleRegulator(sim, sim.moduleName)	# use this file for reg too
    # Only biophys and regulation
    sim.init(biophys, regul, None, None)

    sim.addCell(cellType=0, pos=(0,0,0))

    if sim.is_gui:
        # Add some objects to draw the models
        from CellModeller.GUI import Renderers
        therenderer = Renderers.GLBacteriumRenderer(sim)
        sim.addRenderer(therenderer)

    sim.pickleSteps = 10

def init(cell):
    cell.targetVol = 3.5 + random.uniform(0.0,0.5)
    cell.growthRate = 1.0
    cell.N0 = 0
    cell.p1, cell.p2, cell.p3 = 0,0,0

def repressilator_step(cell, dt):
    t = 0
    #for i in range(10000):
    while t<dt:
        # Solve the quadratic to find free repressor concentrations **n
        p1n = cell.p1**n
        p2n = cell.p2**n
        p3n = cell.p3**n
        p1freen = p1n - Kn - 2*cell.N0 + np.sqrt( Kn*Kn + (2*cell.N0-p1n)**2 + 2*Kn*(2*cell.N0+p1n) )
        p2freen = p2n - Kn - 2*cell.N0 + np.sqrt( Kn*Kn + (2*cell.N0-p2n)**2 + 2*Kn*(2*cell.N0+p2n) )
        p3freen = p3n - Kn - 2*cell.N0 + np.sqrt( Kn*Kn + (2*cell.N0-p3n)**2 + 2*Kn*(2*cell.N0+p3n) )

        # Reaction propensities
        # Production of repressors
        a1_1 = l * cell.N0 * Kn / (Kn + p3freen)
        a1_2 = l * cell.N0 * Kn / (Kn + p1freen)
        a1_3 = l * cell.N0 * Kn / (Kn + p2freen)
        # Degradation of repressors
        a2_1 = cell.p1
        a2_2 = cell.p2
        a2_3 = cell.p3
        # Replication of plasmid
        a3 = _N0_
        # Degradation of plasmid
        a4 = cell.N0

        # Total of propensities
        A = a1_1 + a1_2 + a1_3 + a2_1 + a2_2 + a2_3 + a3 + a4

        # Time to next reaction
        tau = 1/A * np.log(1/np.random.random())
        t += tau

        # Random number to select next reaction
        a_i = np.random.random() * A

        # Random number for protein count update on production
        b = np.random.geometric(1/_b_) # Note _b_ = mean = 1/p

        # Update protein and plasmid numbers according to reaction selected
        if a_i < a1_1:
            # Production of repressor 1
            cell.p1 += b
        elif a_i < a1_1 + a1_2:
            # Production of repressor 2
            cell.p2 += b
        elif a_i < a1_1 + a1_2 + a1_3:
            # Production of repressor 3
            cell.p3 += b
        elif a_i < a1_1 + a1_2 + a1_3 + a2_1:
            # Degradation of repressor 1
            cell.p1 -= 1
        elif a_i < a1_1 + a1_2 + a1_3 + a2_1 + a2_2:
            # Degradation of repressor 2
            cell.p2 -= 1
        elif a_i < a1_1 + a1_2 + a1_3 + a2_1 + a2_2 + a2_3:
            # Degradation of repressor 3
            cell.p3 -= 1
        elif a_i < a1_1 + a1_2 + a1_3 + a2_1 + a2_2 + a2_3 + a3:
            # Replication of plasmid
            cell.N0 += 1
        else:
            # Degradation of plasmid
            cell.N0 -= 1

def update(cells):
    for (id, cell) in cells.items():
        cell.color = [0.1, cell.n_a/20.0, cell.n_b/20.0]
        if cell.volume > cell.targetVol:
            cell.divideFlag = True
        repressilator_step(cell, 0.01*60)

def divide(parent, d1, d2):
    d1.targetVol = 3.5 + random.uniform(0.0,0.5)
    d2.targetVol = 3.5 + random.uniform(0.0,0.5)

    # Binomial separation of plasmids and proteins
    d1.N0 = np.random.binomial(parent.N0, 0.5)
    d2.N0 = parent.N0 - d1.N0
    d1.p1 = np.random.binomial(parent.p1, 0.5)
    d2.p1 = parent.p1 - d1.p1
    d1.p2 = np.random.binomial(parent.p2, 0.5)
    d2.p2 = parent.p2 - d1.p2
    d1.p3 = np.random.binomial(parent.p3, 0.5)
    d2.p3 = parent.p3 - d1.p3
