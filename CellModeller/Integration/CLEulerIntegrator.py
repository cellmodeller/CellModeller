import numpy
import scipy.integrate.odepack
import pyopencl as cl
import pyopencl.array as cl_array
from pyopencl.array import vec
import math

class CLEulerIntegrator:
    #Simple forward Euler integration of species rates
    
    def __init__(self, sim, nSpecies, maxCells, regul=None):
        self.sim = sim
        self.dt = self.sim.dt
        self.regul = regul

        self.cellStates = sim.cellStates
        self.nCells = len(self.cellStates)

        self.nSpecies = nSpecies
        self.maxCells = maxCells

        self.maxSpecDataLen = self.maxCells*nSpecies
        # no need to scale up signal storage
        storageLen = self.maxSpecDataLen

        # These arrays store the level and rate of signals and species
        # in a contiguous form. 
        # To avoid reallocation, create enough space for maxCells
        self.levels = numpy.zeros(storageLen,dtype=numpy.float32)
        self.rates = numpy.zeros(storageLen,dtype=numpy.float32)
        self.makeViews()

        (self.context, self.queue) = self.sim.getOpenCL()
        self.initArrays()
        #self.initKernels()

        # set the species for existing states to views of the levels array
        cs = self.cellStates
        for id,c in list(cs.items()):
            c.species = self.specLevel[c.idx,:]


    def makeViews(self):
        # Level views (references) to the data
        self.specLevel = self.levels[0:self.maxSpecDataLen].reshape(self.maxCells,self.nSpecies)
        # Rate views (references) to the data
        self.specRate = self.rates[0:self.maxSpecDataLen].reshape(self.maxCells,self.nSpecies)

    def addCell(self, cellState):
        idx = cellState.idx
        self.nCells += 1
        cellState.species = self.specLevel[idx,:]
        self.celltype[idx] = numpy.int32(cellState.cellType)

    def divide(self, pState, d1State, d2State):
        # Simulator should have organised indexing:

        # Set up slicing of levels for each daughter and copy parent levels
        d1idx = d1State.idx
        self.nCells += 1
        self.specLevel[d1idx,:] = pState.species
        d1State.species = self.specLevel[d1idx,:]
        self.celltype[d1idx] = d1State.cellType

        d2idx = d2State.idx
        self.nCells += 1
        self.specLevel[d2idx,:] = pState.species
        d1State.species = self.specLevel[d1idx,:]
        self.celltype[d1idx] = d1State.cellType

        d2idx = d2State.idx
        self.nCells += 1
        self.specLevel[d2idx,:] = pState.species
        d2State.species = self.specLevel[d2idx,:]
        self.celltype[d2idx] = d2State.cellType

    def setRegulator(self, regul):
        self.regul = regul
        # Use regulation module to setup kernels
        self.initKernels()

    def initArrays(self):
        self.specLevel_dev = cl_array.zeros(self.queue, (self.maxCells,self.nSpecies), dtype=numpy.float32)
        self.specRate_dev = cl_array.zeros(self.queue, (self.maxCells,self.nSpecies), dtype=numpy.float32)

        self.celltype = numpy.zeros((self.maxCells,), dtype=numpy.int32)
        self.celltype_dev = cl_array.zeros(self.queue, (self.maxCells,),dtype=numpy.int32)
    
        self.effgrow = numpy.zeros((self.maxCells,), dtype=numpy.float32)
        self.effgrow_dev = cl_array.zeros(self.queue, (self.maxCells,), dtype=numpy.float32)
    
        #self.pos_dev = cl_array.zeros(self.queue, (self.maxCells,), dtype=vec.float4)

    def initKernels(self):
        # Get user defined kernel source
        specRateKernel = self.regul.specRateCL()
        from pkg_resources import resource_string
        kernel_src = resource_string(__name__, 'CLEulerIntegrator.cl').decode()
        # substitute user defined kernel code, and number of signals
        kernel_src = kernel_src%(specRateKernel)
        self.program = cl.Program(self.context, kernel_src).build(cache_dir=False)


    def dydt(self):
        self.celltype_dev.set(self.celltype)
        # compute species rates
        self.specLevel_dev.set(self.specLevel)
        self.program.speciesRates(self.queue, (self.nCells,), None,
                                  numpy.int32(self.nSpecies),
                                  self.sim.phys.cell_centers_dev.data,
                                  self.sim.phys.cell_areas_dev.data,
                                  self.sim.phys.cell_vols_dev.data,
                                  self.celltype_dev.data,
                                  self.effgrow_dev.data,
                                  self.specLevel_dev.data,
                                  self.specRate_dev.data).wait()
        self.specRate[:] = self.specRate_dev.get()

    def step(self, dt):
        if dt!=self.dt:
            print("I can only integrate at fixed dt!")
            return

        self.nCells = len(self.cellStates)
        # Check we have enough space allocated
        try:
            s = self.specLevel[self.nCells-1]
        except IndexError:
            # Could resize here, then would have to rebuild views
            print("Number of cells exceeded " \
                    + self.__class__.__name__ \
                    + "::maxCells (" + self.maxCells + ")")

        self.dataLen = self.nCells*self.nSpecies

        self.cellStates = self.sim.cellStates
        cs = self.cellStates
        for id,c in list(cs.items()):
            self.effgrow[c.idx] = numpy.float32(c.effGrowth)
        self.effgrow_dev.set(self.effgrow)

        # growth dilution of species
        self.diluteSpecies()

        self.dydt()
        self.rates[0:self.dataLen] *= self.dt
        self.levels[0:self.dataLen] += self.rates[0:self.dataLen]


# Put the final signal levels into the cell states
#        states = self.cellStates
#        for (id,c) in states.items():
#            if self.signalling:
#                c.signals = self.signalling.signals(c, self.signalLevel)

    def setLevels(self, specLevel):
        self.cellStates = self.sim.cellStates
        self.levels = specLevel
        self.makeViews()
        self.specLevel_dev.set(self.specLevel)
        cs = self.cellStates
        for id,c in list(cs.items()):
            c.species = self.specLevel[c.idx,:]
            self.celltype[c.idx] = numpy.int32(c.cellType)
        self.celltype_dev.set(self.celltype)


    def diluteSpecies(self):
        self.specLevel_dev.set(self.specLevel)
        self.program.diluteSpecs(self.queue, (self.nCells,), None,
                                 numpy.int32(self.nSpecies),
                                 self.sim.phys.cell_old_vols_dev.data,
                                 self.sim.phys.cell_vols_dev.data,
                                 self.specLevel_dev.data).wait()
        self.specLevel[:] = self.specLevel_dev.get()

