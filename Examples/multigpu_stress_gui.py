"""Single-founder colony growth for the standard CellModeller GUI.

Use Load Model and select this file, choose the GPU(s), then click Run.
The setup follows ex1_simpleGrowth2D.py: one cell at the origin, ordinary
2D rod mechanics, and growth/division into a colony. Edit the constants below
for larger runs. Selecting GPUs distributes work; it does not guarantee full
utilization. Use multigpu_stress.py for headless diagnostics and species tests.
"""

import random
from multigpu_stress import capacity_estimate, host_available_bytes

from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium

# Auto estimates a target from host RAM and selected GPUs; an integer overrides it.
max_cells = "auto"
memory_fraction = 0.7
max_contacts = 32
max_sqs = None  # Automatically sized for the selected target.
random_seed = 12345
growth_rate = 2.0
report_steps = 100


def setup(sim):
    global _sim, _rng, _capacity
    _sim = sim
    _rng = random.Random(random_seed)
    random.seed(random_seed)

    _capacity = (capacity_estimate(sim, max_contacts, 0, memory_fraction,
                                   host_available_bytes())
                 if max_cells == "auto" else int(max_cells))
    if _capacity < 1 or max_contacts < 8 or _capacity * max_contacts * 8 >= 2**31:
        raise ValueError('Invalid population/contact capacity')
    grid_capacity = max_sqs if max_sqs is not None else max(192**2, 4 * _capacity)
    print('Colony capacity target: %d cells; grid: %d bins. '
          'This is a memory estimate, not a measured maximum.' % (_capacity, grid_capacity))

    class CapacityPhysics(CLBacterium):
        def update_grid(self):
            super().update_grid()
            if self.n_sqs > self.max_sqs:
                raise MemoryError('Spatial grid exhausted; increase max_sqs in the model')

    biophys = CapacityPhysics(sim, jitter_z=False, max_cells=_capacity,
                             max_contacts=max_contacts, max_sqs=grid_capacity)
    regul = ModuleRegulator(sim)
    sim.init(biophys, regul, None, None)
    sim.addCell(cellType=0, pos=(0, 0, 0))

    if sim.is_gui:
        from CellModeller.GUI import Renderers
        sim.addRenderer(Renderers.GLBacteriumRenderer(sim))

    # Use the normal Save Pickles control when output is needed.
    sim.pickleSteps = 100
    sim.setSaveOutput(False)
    sim.CLWorkStats.clear()


def init(cell):
    cell.targetVol = 2.5 + _rng.uniform(0.0, 0.5)
    cell.growthRate = growth_rate
    cell.color = [0.1, 1.0, 0.3]


def update(cells):
    # Each division adds one cell. Reserve a slot per parent to avoid an
    # overflowing final division wave, then stop growth at the capacity.
    slots = max(0, _capacity - len(cells))
    for cell in cells.values():
        cell.divideFlag = bool(slots and cell.volume > cell.targetVol)
        if cell.divideFlag:
            slots -= 1
        cell.growthRate = growth_rate
    if slots == 0:
        for cell in cells.values():
            cell.growthRate = 0.0
    if _sim.stepNum and _sim.stepNum % report_steps == 0:
        work = _sim.CLWorkStats.get('physics.find_contacts')
        print('Colony: %d / %d cells; cumulative contact work per GPU: %s'
              % (len(cells), _capacity, work if work is not None else 'single device'))


def divide(parent, d1, d2):
    d1.targetVol = 2.5 + _rng.uniform(0.0, 0.5)
    d2.targetVol = 2.5 + _rng.uniform(0.0, 0.5)
