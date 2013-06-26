import numpy
import scipy.integrate.odepack
from scipy.sparse.linalg import LinearOperator
from scipy.ndimage.filters import convolve
from scipy.sparse.linalg import gmres
import pyopencl as cl
import pyopencl.array as cl_array
from pyopencl.array import vec
import math

class CrankNicIntegrator:
    def __init__(self, sim, nSignals, nSpecies, maxCells, sig, greensThreshold=1e-12, regul=None):
        self.sim = sim
        self.dt = self.sim.dt
        self.greensThreshold = greensThreshold
        self.regul = regul
        self.cellStates = sim.cellStates
        self.nCells = len(self.cellStates)

        self.nSpecies = nSpecies
        self.nSignals = nSignals
        self.maxCells = maxCells

        # The signalling model, must be a grid based thing
        self.signalling = sig
        self.gridDim = sig.gridDim

        self.signalDataLen = self.signalling.dataLen()
        self.maxSpecDataLen = self.maxCells*nSpecies
        # no need to scale up signal storage
        storageLen = self.maxSpecDataLen + self.signalDataLen

        # These arrays store the level and rate of signals and species
        # in a contiguous form. The first part is the signals,
        # then the cell species
        # To avoid reallocation, create enough space for maxCells

        self.levels = numpy.zeros(storageLen)
        self.rates = numpy.zeros(storageLen)
        self.makeViews()
        # Set initial distribution of signals
        if self.signalling.initLevels:
            for s in range(self.nSignals):
                grid = self.signalLevel.reshape(self.gridDim)
                grid[s,:] = self.signalling.initLevels[s]

        self.computeGreensFunc()

        # Initialise map of cell ids to index in arrays
        # set the species for existing states to views of the levels array
        cs = self.cellStates
        for c in cs.items():
            c.species = self.specLevels[c.idx,:]

    def makeViews(self):
        # Level views (references) to the data
        self.signalLevel = self.levels[0:self.signalDataLen]
        self.specLevel = self.levels[self.signalDataLen:self.signalDataLen+self.maxSpecDataLen].reshape(self.maxCells,self.nSpecies)
        # Rate views (references) to the data
        self.signalRate = self.rates[0:self.signalDataLen]
        self.specRate = self.rates[self.signalDataLen:self.signalDataLen+self.maxSpecDataLen].reshape(self.maxCells,self.nSpecies)

    def CNOperator(self, v):
        # Transport operator
        self.signalling.transportRates(self.signalRate, v)
        # Return (I-hT/2)v, where T is transport operator, h=dt
        return v - 0.5*self.dt*self.signalRate

    def computeGreensFunc(self):
        L = LinearOperator((self.signalDataLen,self.signalDataLen), matvec=self.CNOperator, dtype=numpy.float32)
        rhs = numpy.zeros(self.gridDim, dtype=numpy.float32)
        idx = ( math.floor(self.gridDim[1]*0.5), math.floor(self.gridDim[2]*0.5), math.floor(self.gridDim[3]*0.5) )
        for s in xrange(self.nSignals):
            rhs[(s,)+idx] = 1.0 # ~delta function in each signal
        (self.greensFunc, info) = gmres(L,rhs.reshape(self.signalDataLen)) # Solve impulse response = greens func
        # Take only bounding box of region where G > threshold
        self.greensFunc.shape = self.gridDim
        inds = numpy.transpose(numpy.nonzero(self.greensFunc.reshape(self.gridDim)>self.greensThreshold))
        self.greensFunc = self.greensFunc[:, min(inds[:,1]):max(inds[:,1])+1, \
                                            min(inds[:,2]):max(inds[:,2])+1, \
                                            min(inds[:,3]):max(inds[:,3])+1]
        print "Truncated Green's function size is " + str(self.greensFunc.shape)


    def addCell(self, cellState):
        idx = cellState.idx
        self.nCells += 1
        cellState.species = self.specLevel[idx,:]

    def divide(self, pState, d1State, d2State):
        # Simulator should have organised indexing:

        # Set up slicing of levels for each daughter and copy parent levels
        d1idx = d1State.idx
        self.nCells += 1
        self.specLevel[d1idx,:] = pState.species
        d1State.species = self.specLevel[d1idx,:]

        d2idx = d2State.idx
        self.nCells += 1
        self.specLevel[d2idx,:] = pState.species
        d2State.species = self.specLevel[d2idx,:]

    def setSignalling(self, sig):
        self.sig = sig

    def setRegulator(self, regul):
        self.regul = regul

    def dydt(self):
        # compute cell species production rates into rates array

        # Loop over cells to get rates
        states = self.cellStates
        for (id,c) in states.items():
            idx = c.idx
            cellSignals = self.signalling.signals(c, self.signalLevel)
            self.specRate[idx,:] = self.regul.speciesRates(c, self.specLevel[idx,:], cellSignals)
            cellRates = self.regul.signalRates(c, self.specLevel[idx,:], 
            cellSignals)
            self.signalling.cellProdRates(self.signalRate, c, cellRates)

    def step(self, dt):
        if dt!=self.dt:
            print "I can only integrate at fixed dt!"
            return

        self.nCells = len(self.cellStates)
        # Check we have enough space allocated
        try:
            s = self.specLevel[self.nCells-1]
        except IndexError:
            # Could resize here, then would have to rebuild views
            print "Number of cells exceeded " \
                    + self.__class__.__name__ \
                    + "::maxCells (" + self.maxCells + ")"

        self.dataLen = self.signalDataLen + self.nCells*self.nSpecies

        # Do u += h(T(u_t)/2 + hf(u_t)) where T=transport operator, f(u_t) is 
        # our regulation function dydt
        self.signalling.transportRates(self.signalRate, self.signalLevel)
        self.signalRate *= 0.5
        self.dydt()
        self.rates[0:self.dataLen] *= self.dt
        self.levels[0:self.dataLen] += self.rates[0:self.dataLen]

        # Convolve (I+hT/2)u_t + f(u_t) with the Greens func to get u_{t+1}
        sigLvl = self.signalLevel.reshape(self.gridDim)
        convolve(sigLvl, self.greensFunc, mode='nearest')

        # Put the final signal levels into the cell states
        states = self.cellStates
        for (id,c) in states.items():
            if self.signalling:
                c.signals = self.signalling.signals(c, self.signalLevel)
