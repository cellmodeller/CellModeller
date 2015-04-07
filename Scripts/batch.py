import os
import sys
import subprocess
import string
import shutil

from CellModeller.Simulator import Simulator

max_cells = 100000
cell_buffer = 256

def simulate(modfilename, steps=50):
    (path,name) = os.path.split(modfilename)
    modname = str(name).split('.')[0]
    sys.path.append(path)
    sim = Simulator(modname, 0.25, None)
    while len(sim.cellStates) < max_cells-cell_buffer:
        sim.step()

simulate(sys.argv[1])
