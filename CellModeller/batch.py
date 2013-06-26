import os
import sys
import subprocess
import string
import shutil

from CellModeller.Simulator import Simulator

sys.path.append('../../Models')

max_cells = 30000
cell_buffer = 256

def simulate(mod_name, steps=50):
    sim = Simulator(mod_name, 0.25)
    while len(sim.cellStates) < max_cells-cell_buffer:
        sim.step()

simulate(sys.argv[1])
