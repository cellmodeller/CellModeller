import os
import sys
import subprocess
import string
import shutil

from CellModeller.Simulator import Simulator

sys.path.append('../../Models')

max_cells = 100000
cell_buffer = 256

gammas = [2,5,10,15,30,50,100,500]

def simulate(mod_name, ingam, fname):
    sim = Simulator(mod_name, 0.25, ingam, 25, fname)
    while len(sim.cellStates) < max_cells-cell_buffer:
        sim.step()

#here we are varying gamma from 1 to 50, every 2
for i in range(len(gammas)):
    simulate(sys.argv[1],gammas[i],(str(sys.argv[1])+'_gamma='+str(gammas[i]))
