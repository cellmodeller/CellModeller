import sys
import math
import numpy
import pyopencl as cl
import pyopencl.array as cl_array
from pyopencl.array import vec

class CLFixedPosition:
    def __init__(self, sim, max_cells=100):
        self.simulator = sim
        self.max_cells = max_cells
        self.init_cl()
        self.init_data()
        self.n_cells = 0

    def setRegulator(self, reg):
        self.reg = reg
                    
    def init_cl(self):
        if self.simulator:
            (self.context, self.queue) = self.simulator.getOpenCL()

    def init_data(self):
        """Set up the data OpenCL will store on the device."""
        # cell data
        cell_geom = (self.max_cells,)
        self.cell_centers = numpy.zeros(cell_geom, vec.float4)
        self.cell_centers_dev = cl_array.zeros(self.queue, cell_geom, vec.float4)

        # cell geometry
        self.cell_areas_dev = cl_array.zeros(self.queue, cell_geom, numpy.float32)
        self.cell_areas = numpy.zeros(cell_geom, numpy.float32)
        self.cell_vols_dev = cl_array.zeros(self.queue, cell_geom, numpy.float32)
        self.cell_vols = numpy.zeros(cell_geom, numpy.float32)
        self.cell_old_vols_dev = cl_array.zeros(self.queue, cell_geom, numpy.float32)
        self.cell_old_vols = numpy.zeros(cell_geom, numpy.float32)

        self.cell_growth_rates = numpy.zeros(cell_geom, numpy.float32)
        self.cell_growth_rates_dev = cl_array.zeros(self.queue, cell_geom, numpy.float32)

    def addCell(self, cellState, pos=(0,0,0), vol=1, **kwargs):
        i = cellState.idx
        self.n_cells += 1
        cid = cellState.id
        self.cell_centers[i] = tuple(pos+(0,))
        self.cell_vols[i] = vol
        self.cell_old_vols[i] = vol
        self.initCellState(cellState)
        self.set_cells()

    def initCellState(self, state):
        cid = state.id
        i = state.idx
        state.pos = [self.cell_centers[i][j] for j in range(3)]
        state.volume = self.cell_vols[i]
        state.startVol = state.volume

    def updateCellState(self, state):
        cid = state.id
        i = state.idx
        state.volume = self.cell_vols[i] 
        self.cell_growth_rates[i] = state.growthRate*state.volume

    def get_cells(self):
        """Copy cell centers, dirs, lens, and rads from the device."""
        self.cell_centers[0:self.n_cells] = self.cell_centers_dev[0:self.n_cells].get()
        self.cell_vols[0:self.n_cells] = self.cell_vols_dev[0:self.n_cells].get()

    def set_cells(self):
        """Copy cell centers, dirs, lens, and rads to the device from local."""
        self.cell_centers_dev[0:self.n_cells].set(self.cell_centers[0:self.n_cells])
        self.cell_vols_dev[0:self.n_cells].set(self.cell_vols[0:self.n_cells])
        self.cell_growth_rates_dev[0:self.n_cells].set(self.cell_growth_rates[0:self.n_cells])

    def step(self, dt):
        self.set_cells()

        # No change in positions, just growth increasing volume
        self.cell_old_vols_dev[0:self.n_cells] = self.cell_vols_dev[0:self.n_cells]
        self.cell_vols_dev += self.cell_growth_rates_dev * dt

        if self.simulator:
            self.get_cells()
            for state in self.simulator.cellStates.values():
                self.updateCellState(state)

        return True

