import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.GeneralModels.CLSPP import CLSPP
from CellModeller.GUI import Renderers
import numpy as np
import math

#Import Euler integrator for solving ODE system of chemical species inside the cells
from CellModeller.Integration.CLEulerIntegrator import CLEulerIntegrator

n_cells = 1000
side = 40
Wc = 0.5
psi = 0

Fm = 0.5
gamma_s = 1
D = 1. 
fcil = 2 * D * psi

def setup(sim):
    sim.dt = 0.1
    # Set biophysics, signalling, and regulation models
    
    #max_sqr
    #biophys = CLBacterium(sim, max_cells=max_cells, jitter_z=False, max_sqs=256**2)
    biophys = CLSPP(sim, 
            max_cells=n_cells, 
            gamma_s=gamma_s, 
            Fm=Fm,
            Wc=Wc,
            fcil=fcil,
            D=D,
            max_planes=6, 
            grid_spacing=2,
            cgs_tol=1e-4,
            max_substeps=1,
            spherical=False)

    # use this file for reg too
    regul = ModuleRegulator(sim)
    # Only biophys and regulation
    sim.init(biophys, regul, None, None)

    # Specify the initial cell and its location in the simulation
    # Grid
    '''
    for i in range(-20,20):
        for j in range(-20,20):
            p = (i*2, j*2, 0) # np.random.uniform(-1,1)
            cell_dir = np.random.uniform(-1,1, size=(3,))
            cell_dir[2] = 0
            cell_dir /= np.sqrt(np.sum(cell_dir*cell_dir))
            sim.addCell(cellType=0, dir=tuple(cell_dir), pos=tuple(p))
    '''
    # Random
    for i in range(n_cells):
        p = np.random.uniform(-1,1, size=(3,)) * side
        p[2] = 0
        cell_dir = np.random.uniform(-1,1, size=(3,))
        cell_dir[2] = 0
        cell_dir /= np.sqrt(np.sum(cell_dir*cell_dir))
        sim.addCell(cellType=0, dir=tuple(cell_dir), pos=tuple(p))

    # Box
    biophys.addPlane((0,-side,0),(0,1,0), 20.)
    biophys.addPlane((0,side,0),(0,-1,0), 20.)
    biophys.addPlane((-side,0,0),(1,0,0), 20.)
    biophys.addPlane((side,0,0),(-1,0,0), 20.)

    # Add some objects to draw the models
    therenderer = Renderers.GLSphereRenderer(sim, draw_axis=False, draw_nbr_dir=False)
    sim.addRenderer(therenderer)

    sim.pickleSteps = 10

def init(cell):
    cell.color = [1,1,1]

def update(cells):
    pass

def divide(parent, d1, d2):
    pass

