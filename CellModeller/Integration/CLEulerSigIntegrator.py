import numpy
import scipy.integrate.odepack
from scipy.sparse.linalg import LinearOperator
from scipy.ndimage.filters import convolve
from scipy.sparse.linalg import gmres
import pyopencl as cl
import pyopencl.array as cl_array
from pyopencl.array import vec
import math
from functools import reduce


def unique_stable(ar, return_index=False, return_inverse=False):
    """
    Find the unique elements of an array.

    Copied from numpy's unique, but uses stable sorts.
    """
    import numpy as np
    try:
        ar = ar.flatten()
    except AttributeError:
        if not return_inverse and not return_index:
            items = sorted(set(ar)) # sorted is stable
            return np.asarray(items)
        else:
            ar = np.asanyarray(ar).flatten()

    if ar.size == 0:
        if return_inverse and return_index:
            return ar, np.empty(0, np.bool), np.empty(0, np.bool)
        elif return_inverse or return_index:
            return ar, np.empty(0, np.bool)
        else:
            return ar

    if return_inverse or return_index:
        perm = ar.argsort(kind='mergesort') # using mergesort for stability
        aux = ar[perm]
        flag = np.concatenate(([True], aux[1:] != aux[:-1]))
        if return_inverse:
            iflag = np.cumsum(flag) - 1
            iperm = perm.argsort(kind='mergesort') # using mergesort for stability
            if return_index:
                return aux[flag], perm[flag], iflag[iperm]
            else:
                return aux[flag], iflag[iperm]
        else:
            return aux[flag], perm[flag]

    else:
        ar.sort(kind='mergesort')
        flag = np.concatenate(([True], ar[1:] != ar[:-1]))
        return ar[flag]




class CLEulerSigIntegrator:
    def __init__(self, sim, nSignals, nSpecies, maxCells, sig, deg_rate=0, regul=None, boundcond='constant'):
        self.sim = sim
        self.dt = self.sim.dt
        self.regul = regul
        self.boundcond = boundcond


        self.nSpecies = nSpecies
        self.nSignals = nSignals
        self.maxCells = maxCells

        self.deg_rate = deg_rate

        # The signalling model, must be a grid based thing
        self.signalling = sig
        self.gridDim = sig.gridDim
        self.gridTotalSize = reduce(lambda x, y: x * y, self.gridDim[1:4])

        self.signalDataLen = self.signalling.dataLen()
        self.maxSpecDataLen = self.maxCells*nSpecies
        # no need to scale up signal storage
        storageLen = self.maxSpecDataLen + self.signalDataLen

        # These arrays store the level and rate of signals and species
        # in a contiguous form. The first part is the signals,
        # then the cell species
        # To avoid reallocation, create enough space for maxCells

        self.levels = numpy.zeros(storageLen,dtype=numpy.float32)
        self.rates = numpy.zeros(storageLen,dtype=numpy.float32)
        self.makeViews()

        # Set initial distribution of signals
        if self.signalling.initLevels:
            for s in range(self.nSignals):
                grid = self.signalLevel.reshape(self.gridDim)
                grid[s,:] = self.signalling.initLevels[s]

        (self.context, self.queue) = self.sim.getOpenCL()
        self.initArrays()
        #self.initKernels()

        self.setCellStates(sim.cellStates)
        

    def saveData(self, data):
        integ_data = {
                'specData': self.levels,
                'sigGrid': self.signalLevel,
                'sigData': self.cellSigLevels,
                }
        data.update(integ_data)
        return data

    def loadData(self, data):
        self.levels = data['specData']
        self.signalLevel = data['sigGrid']
        self.cellSigLevels = data['sigData']
        self.setCellStates(self.sim.cellStates)

    '''
    def setLevels(self, SSLevel, cellSigData):
        self.cellStates = self.sim.cellStates
        self.levels = SSLevel
        self.makeViews()
        self.cellSigLevels = cellSigData
        self.signalLevel_dev.set(self.signalLevel)
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
    '''

    def setCellStates(self, cs):
        # set the species for existing states to views of the levels array
        self.cellStates = cs
        self.nCells = len(self.cellStates)
        self.cellSigLevels[:,:] = 0
        for id,c in list(cs.items()):
            c.species = self.specLevel[c.idx,:]
            c.signals = self.cellSigLevels[c.idx,:]
            c.gradient = self.cellSigGradients[c.idx,:]
            self.celltype[c.idx] = numpy.int32(c.cellType)


    def makeViews(self):
        # Level views (references) to the data
        self.signalLevel = self.levels[0:self.signalDataLen]
        self.specLevel = self.levels[self.signalDataLen:self.signalDataLen+self.maxSpecDataLen].reshape(self.maxCells,self.nSpecies)
        # Rate views (references) to the data
        self.signalRate = self.rates[0:self.signalDataLen]
        self.specRate = self.rates[self.signalDataLen:self.signalDataLen+self.maxSpecDataLen].reshape(self.maxCells,self.nSpecies)

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
        self.gridIdxs = numpy.zeros((self.maxCells,8),dtype=numpy.int32)
        self.gridIdxs_dev = cl_array.zeros(self.queue, (self.maxCells,8),dtype=numpy.int32)
        self.triWts = numpy.zeros((self.maxCells,8),dtype=numpy.float32)
        self.triWts_dev = cl_array.zeros(self.queue, (self.maxCells,8),dtype=numpy.float32)
        self.cellSigRates = numpy.zeros((self.maxCells,8,self.nSignals),dtype=numpy.float32)
        self.cellSigRates_dev = cl_array.zeros(self.queue, (self.maxCells,8,self.nSignals),dtype=numpy.float32)
        self.cellSigLevels = numpy.zeros((self.maxCells,self.nSignals),dtype=numpy.float32)
        self.cellSigLevels_dev = cl_array.zeros(self.queue, (self.maxCells,self.nSignals),dtype=numpy.float32)
        self.cellSigGradients = numpy.zeros((self.maxCells,self.nSignals),dtype=vec.float4)
        self.cellSigGradients_dev = cl_array.zeros(self.queue, (self.maxCells,self.nSignals),dtype=vec.float4)
        self.signalLevel_dev = cl_array.zeros(self.queue, self.gridDim, dtype=numpy.float32)
        self.signalGradient_x_dev = cl_array.zeros(self.queue, self.gridDim, dtype=numpy.float32)
        self.signalGradient_y_dev = cl_array.zeros(self.queue, self.gridDim, dtype=numpy.float32)
        self.signalGradient_z_dev = cl_array.zeros(self.queue, self.gridDim, dtype=numpy.float32)
        self.signalGradient_x = numpy.zeros(self.gridDim, dtype=numpy.float32) 
        self.signalGradient_y = numpy.zeros(self.gridDim, dtype=numpy.float32)
        self.signalGradient_z = numpy.zeros(self.gridDim, dtype=numpy.float32)
        self.specLevel_dev = cl_array.zeros(self.queue, (self.maxCells,self.nSpecies), dtype=numpy.float32)
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
        kernel_src = resource_string(__name__, 'CLEulerSigIntegrator.cl').decode()
        # substitute user defined kernel code, and number of signals
        kernel_src = kernel_src % {'sigKernel': sigRateKernel,
                                   'specKernel': specRateKernel,
                                   'nSignals': self.nSignals}
        self.program = cl.Program(self.context, kernel_src).build(cache_dir=False)


    def dydt(self):
        # get cell grid idxs and weights
        self.program.gridCells(self.queue, (self.nCells,), None,
                numpy.float32(self.signalling.gridOrig[0]),
                numpy.float32(self.signalling.gridOrig[1]),
                numpy.float32(self.signalling.gridOrig[2]),
                numpy.float32(self.signalling.gridSize[0]),
                numpy.float32(self.signalling.gridSize[1]),
                numpy.float32(self.signalling.gridSize[2]),
                numpy.int32(self.signalling.gridDim[1]),
                numpy.int32(self.signalling.gridDim[2]),
                numpy.int32(self.signalling.gridDim[3]),
                self.sim.phys.cell_centers_dev.data,
                self.triWts_dev.data,
                self.gridIdxs_dev.data).wait()
        self.gridIdxs[:] = self.gridIdxs_dev.get()

        # put local cell signal levels in array
        self.signalLevel_dev.set(self.signalLevel)
        self.program.setCellSignals(self.queue, (self.nCells,), None,
                numpy.int32(self.nSignals),
                numpy.int32(self.gridTotalSize),
                numpy.int32(self.signalling.gridDim[1]),
                numpy.int32(self.signalling.gridDim[2]),
                numpy.int32(self.signalling.gridDim[3]),
                self.gridIdxs_dev.data,
                self.triWts_dev.data,
                self.signalLevel_dev.data,
                self.signalGradient_x_dev.data,
                self.signalGradient_y_dev.data,
                self.signalGradient_z_dev.data,
                self.cellSigLevels_dev.data,
                self.cellSigGradients_dev.data
                ).wait()

        self.celltype_dev.set(self.celltype)
        # compute species rates
        self.specLevel_dev.set(self.specLevel)
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
        self.specRate[:] = self.specRate_dev.get()

        # compute signal rates, weighted for grid nodes
        self.program.signalRates(self.queue, (self.nCells,), None,
                                 numpy.int32(self.nSignals),
                                 numpy.int32(self.nSpecies),
                                 numpy.float32(self.sim.sig.dV),
                                 self.sim.phys.cell_areas_dev.data,
                                 self.sim.phys.cell_vols_dev.data,
                                 self.celltype_dev.data,
                                 self.specLevel_dev.data,
                                 self.cellSigLevels_dev.data,
                                 self.triWts_dev.data,
                                 self.cellSigRates_dev.data).wait()
        self.cellSigRates[:] = self.cellSigRates_dev.get()

        # Put cell signal production into diffusion grid:
        #  - Using a convoluted way to reduce by key to get cell prod rates
        #  into grid
        gridIdxs = self.gridIdxs[0:self.nCells,:].reshape(self.nCells*8)
        order = numpy.argsort(gridIdxs)
        gridIdxs = gridIdxs[order]
        cellSigRates = self.cellSigRates[0:self.nCells,:,:].reshape((self.nCells*8,self.nSignals))
        cellSigRates = cellSigRates[order]
        cellSigRates.cumsum(axis=0, out=cellSigRates)
        (u,index) = unique_stable(gridIdxs[::-1],True) # first occurance of each key in reverse keys
        index = len(gridIdxs)-1-index # last occurance in forward keys
        cellSigRates = cellSigRates[index]
        idxs = gridIdxs[index]
        cellSigRates[1:] = cellSigRates[1:] - cellSigRates[:-1] # difference in cumsums is sum for each index


        signalRate = self.signalRate.reshape((self.nSignals,self.gridDim[1]*self.gridDim[2]*self.gridDim[3]))
        signalRate[:,idxs] += cellSigRates.transpose() # add into diffusion grid

        sigrate = signalRate.reshape(self.gridDim)
        x = int(self.gridDim[1]/2)
        y = int(self.gridDim[2]/2)
        z = int(self.gridDim[3]/2)
        sigrate[:,x,y,z] += 100.


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

        self.dataLen = self.signalDataLen + self.nCells*self.nSpecies

        # growth dilution of species
        self.diluteSpecies()

        # Laplacian
        self.signalling.transportRates(self.signalRate, self.signalLevel, self.boundcond)
        # Degradation
        self.signalRate -= self.signalLevel * self.deg_rate
        # Gradient
        self.signalling.gradient(self.signalGradient_x, self.signalGradient_y, self.signalGradient_z, self.signalLevel, self.boundcond)
        self.dydt()
        self.rates[0:self.dataLen] *= self.dt
        self.levels[0:self.dataLen] += self.rates[0:self.dataLen]

        # put local cell signal levels in array
        self.signalLevel_dev.set(self.signalLevel)
        self.signalGradient_x_dev.set(self.signalGradient_x)
        self.signalGradient_y_dev.set(self.signalGradient_y)
        self.signalGradient_z_dev.set(self.signalGradient_z)
        self.program.setCellSignals(self.queue, (self.nCells,), None,
                numpy.int32(self.nSignals),
                numpy.int32(self.gridTotalSize),
                numpy.int32(self.signalling.gridDim[1]),
                numpy.int32(self.signalling.gridDim[2]),
                numpy.int32(self.signalling.gridDim[3]),
                self.gridIdxs_dev.data,
                self.triWts_dev.data,
                self.signalLevel_dev.data,
                self.signalGradient_x_dev.data,
                self.signalGradient_y_dev.data,
                self.signalGradient_z_dev.data,
                self.cellSigLevels_dev.data,
                self.cellSigGradients_dev.data).wait()
        self.cellSigLevels[:] = self.cellSigLevels_dev.get()
        self.cellSigGradients[:] = self.cellSigGradients_dev.get()

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
        self.levels = SSLevel
        self.makeViews()
        self.cellSigLevels = cellSigData
        self.signalLevel_dev.set(self.signalLevel)
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


