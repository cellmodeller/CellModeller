import numpy
import scipy.integrate.odepack
from scipy.sparse.linalg import LinearOperator
from scipy.ndimage.filters import convolve
from scipy.sparse.linalg import gmres
import pyopencl as cl
import pyopencl.array as cl_array
from pyopencl.array import vec
import math
from fenics import *

class CLEulerFEMIntegrator:
    def __init__(self, sim, nSignals, nSpecies, maxCells, sig, deg_rate=0, regul=None):
        self.sim = sim
        self.dt = self.sim.dt
        self.regul = regul

        self.nSpecies = nSpecies
        self.nSignals = nSignals
        self.maxCells = maxCells

        self.deg_rate = deg_rate

        (self.context, self.queue) = self.sim.getOpenCL()
        self.initArrays()

        self.signalling = sig
        self.levels = None
        self.setCellStates(sim.cellStates)

    def saveData(self, data):
        integ_data = {
                'specData': self.specLevel,
                'sigData': self.cellSigLevels,
                }
        data.update(integ_data)
        return data 

    def loadData(self, data):
        self.specLevel = data['specLevel']
        self.cellSigLevels = data['sigData']
        self.specLevel_dev.set(self.specLevel)
        self.cellSigLevels_dev.set(self.cellSigLevels)
        self.setCellStates(self.sim.cellStates)

    def setCellStates(self, cs):
        # set the species for existing states to views of the levels array
        self.cellStates = cs
        self.nCells = len(self.cellStates)
        for id,c in list(cs.items()):
            cellState.species = self.specLevel[c.idx,:]
            cellState.signals = self.cellSigLevels[c.idx,:]
            c.gradient = self.cellSigGradients[c.idx,:]
            self.celltype[c.idx] = numpy.int32(c.cellType)

    def addCell(self, cellState):
        idx = cellState.idx
        self.nCells += 1
        cellState.species = self.specLevel[idx,:]
        cellState.signals = self.cellSigLevels[idx,:]
        cellState.gradient = self.cellSigGradients[idx,:]
        self.celltype[idx] = numpy.int32(cellState.cellType)

    def divide(self, pState, d1State, d2State):
        # Simulator should have organised indexing:

        # Set up slicing of levels for each daughter and copy parent levels
        d1idx = d1State.idx
        self.nCells += 1
        self.specLevel[d1idx,:] = pState.species
        self.cellSigLevels[d1idx,:] = pState.signals
        d1State.species = self.specLevel[d1idx,:]
        d1State.signals = self.cellSigLevels[d1idx,:]
        self.celltype[d1idx] = d1State.cellType

        d2idx = d2State.idx
        self.nCells += 1
        self.specLevel[d2idx,:] = pState.species
        self.cellSigLevels[d2idx,:] = pState.signals
        d1State.species = self.specLevel[d1idx,:]
        d1State.signals = self.cellSigLevels[d1idx,:]
        self.celltype[d1idx] = d1State.cellType

        d2idx = d2State.idx
        self.nCells += 1
        self.specLevel[d2idx,:] = pState.species
        d2State.species = self.specLevel[d2idx,:]
        d2State.signals = self.cellSigLevels[d2idx,:]
        self.celltype[d2idx] = d2State.cellType

    def setSignalling(self, sig):
        self.signalling = sig

    def setRegulator(self, regul):
        self.regul = regul
        # Use regulation module to setup kernels
        self.initKernels()

    def initArrays(self):
        self.cellSigRates = numpy.zeros((self.maxCells,self.nSignals),dtype=numpy.float32)
        self.cellSigRates_dev = cl_array.zeros(self.queue, (self.maxCells,self.nSignals),dtype=numpy.float32)
        self.cellSigLevels = numpy.zeros((self.maxCells,self.nSignals),dtype=numpy.float32)
        self.cellSigLevels_dev = cl_array.zeros(self.queue, (self.maxCells,self.nSignals),dtype=numpy.float32)
        self.cellSigGradients = numpy.zeros((self.maxCells,self.nSignals),dtype=vec.float4)
        self.cellSigGradients_dev = cl_array.zeros(self.queue, (self.maxCells,self.nSignals),dtype=vec.float4)
        self.gx_dev = cl_array.zeros(self.queue, (self.maxCells,self.nSignals),dtype=numpy.float32)
        self.gy_dev = cl_array.zeros(self.queue, (self.maxCells,self.nSignals),dtype=numpy.float32)
        self.gz_dev = cl_array.zeros(self.queue, (self.maxCells,self.nSignals),dtype=numpy.float32)

        self.specLevel = numpy.zeros((self.maxCells,self.nSpecies), dtype=numpy.float32)
        self.specLevel_dev = cl_array.zeros(self.queue, (self.maxCells,self.nSpecies), dtype=numpy.float32)
        self.specRate = numpy.zeros((self.maxCells,self.nSpecies), dtype=numpy.float32)
        self.specRate_dev = cl_array.zeros(self.queue, (self.maxCells,self.nSpecies), dtype=numpy.float32)

        self.celltype = numpy.zeros((self.maxCells,),dtype=numpy.int32)
        self.celltype_dev = cl_array.zeros(self.queue, (self.maxCells,),dtype=numpy.int32)
        #self.pos_dev = cl_array.zeros(self.queue, (self.maxCells,), dtype=vec.float4)

    def initKernels(self):
        # Get user defined kernel source
        specRateKernel = self.regul.specRateCL()
        sigRateKernel = self.regul.sigRateCL()
        #kernel_src = open('CellModeller/Integration/CLCrankNicIntegrator.cl', 'r').read()
        from pkg_resources import resource_string
        kernel_src = resource_string(__name__, 'CLEulerFEMIntegrator.cl').decode()
        # substitute user defined kernel code, and number of signals
        kernel_src = kernel_src % {'sigKernel': sigRateKernel,
                                   'specKernel': specRateKernel,
                                   'nSignals': self.nSignals}
        self.program = cl.Program(self.context, kernel_src).build(cache_dir=False)


    def dydt(self):

        # put local cell signal levels in array
        self.signalling.signals(self.cellSigLevels_dev)
        self.signalling.gradient(self.gx_dev, self.gy_dev, self.gz_dev)
        self.program.combine_grad_components(self.queue, (self.nCells,), None,
                                        self.gx_dev.data,
                                        self.gy_dev.data,
                                        self.gz_dev.data,
         				self.cellSigGradients_dev.data).wait()
        self.celltype_dev.set(self.celltype)
        # compute species rates
        self.program.speciesRates(self.queue, (self.nCells,), None,
                                  numpy.int32(self.nSignals),
                                  numpy.int32(self.nSpecies),
                                  numpy.float32(self.sim.sig.dV),
                                  self.sim.phys.cell_areas_dev.data,
                                  self.sim.phys.cell_vols_dev.data,
                                  self.celltype_dev.data,
                                  self.specLevel_dev.data,
                                  self.cellSigLevels_dev.data,
                                  self.specRate_dev.data).wait()


        self.program.signalRates(self.queue, (self.nCells,), None,
                                 numpy.int32(self.nSignals),
                                 numpy.int32(self.nSpecies),
                                 numpy.float32(self.sim.sig.dV),
                                 self.sim.phys.cell_areas_dev.data,
                                 self.sim.phys.cell_vols_dev.data,
                                 self.celltype_dev.data,
                                 self.specLevel_dev.data,
                                 self.cellSigLevels_dev.data,
                                 self.cellSigRates_dev.data).wait()

        # Get the data back from device
        self.specLevel[:] = self.specLevel_dev.get()
        self.cellSigLevels[:] = self.cellSigLevels_dev.get()
        self.cellSigRates[:] = self.cellSigRates_dev.get()
        self.cellSigGradients[:] = self.cellSigGradients_dev.get()

        # Add point source for each cell 
        for id,c in self.cellStates.items():
            self.signalling.add_point_source(c.pos, self.cellSigRates[c.idx])
        #self.signalling.add_point_source([40,0,0], 100.)

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

        # growth dilution of species
        self.diluteSpecies()

        self.dydt()


# Put the final signal levels into the cell states
#        states = self.cellStates
#        for (id,c) in states.items():
#            if self.signalling:
#                c.signals = self.signalling.signals(c, self.signalLevel)

        # Update cellType array
        for (id,c) in list(self.cellStates.items()):
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

    def setLevels(self, SSLevel, cellSigData):
        self.cellStates = self.sim.cellStates
        self.u = u
        self.cellSigLevels = cellSigData
        self.specLevel_dev.set(self.specLevel)
        self.cellSigLevels_dev.set(self.cellSigLevels)
        self.cellSigGradients_dev.set(self.cellSigGradients)
        cs = self.cellStates
        for id,c in list(cs.items()): #make sure everything is correct here
            c.gradient = self.cellSigGradients[c.idx,:]
            c.species = self.specLevel[c.idx,:]
            c.signals = self.cellSigLevels[c.idx,:]
            self.celltype[c.idx] = numpy.int32(c.cellType)
        self.celltype_dev.set(self.celltype)



