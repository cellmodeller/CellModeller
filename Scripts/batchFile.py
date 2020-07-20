import os
import sys
import subprocess
import string
import shutil
import time
from exceptions import OSError

from CellModeller.Simulator import Simulator

sys.path.append('./')

print((os.getcwd()))

import CellModeller.AdaptiveSimulator


sys.path.append('Models')

cells = 300000
cell_buffer = 256
pickleSteps = 50

try:
    mod_name = sys.argv[1]
except:
    mod_name = 'DaninoOscillator'

try:
    pickleDir = sys.argv[2]
except:
    pickleDir = "/scratch/mg542/batch"

try:
    pickleSteps = int(sys.argv[3])
except:
    'setting pickleSteps to %i' % pickleSteps

try:
    os.mkdir(pickleDir)
except OSError:
    print((pickleDir, 'exists'))

pickleSetDir = os.path.join(pickleDir, mod_name)

try:
    os.mkdir(pickleSetDir)
except:
    print((pickleSetDir, 'exists'))

folderName = time.strftime('%Y%m%d-%H%M%S', time.localtime())
pickleFileRoot = os.path.join(pickleSetDir, folderName)

max_cells = 300000
cell_buffer = 256


def simulate(mod_name, steps=50):
    print('simulate')
    sim = Simulator(mod_name, 0.25, None, pickleSteps=50, pickleFileRoot=pickleFileRoot)
    print('start')

    sim.phys.set_cells()
    sim.phys.calc_cell_geom()

    while len(sim.cellStates) < max_cells-cell_buffer:
        sim.step()

simulate(mod_name)
