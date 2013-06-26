from math import floor
import numpy
from scipy.ndimage import laplace
from scipy.ndimage.filters import convolve


class GridDiffusion:
    """ GridDiffusion: 
    Signal levels/rates on the grid are stored by the Integrator, and this class writes
into that array data.
    1. rates() computes the \Del^2 transport operator for grid based diffusion.
    2. signals() returns the local signal level for a given cell, and given signal 
level array.

    Other transport models can be implemented by replacing the rates() function with
something else - e.g. a bulk flow term.
"""
    def __init__(self, sim, nSignals, gridDim, gridSize, gridOrig, D, adv=None, initLevels=None, regul=None):
        # Data organisation
        self.gridDim = (nSignals,) + gridDim
        self.gridDataLen = reduce(lambda x, y: x * y, self.gridDim)

        self.gridSize = gridSize
        self.gridOrig = numpy.array(gridOrig)
        assert gridSize[0]==gridSize[1]==gridSize[2]

        # Scale diffusion/advection rates by gridSize
        h = gridSize[0]
        self.diffRates = [d/(h**2) for d in D]
        self.advRates = [a/h for a in adv] if adv else None
        self.nSignals = nSignals

        # Initial signal level - homogeneous distribution
        self.initLevels = initLevels

        # This is the volume of the grid element. It should really be used
        # to calc local concentration to a cell. Need to go through and make
        # everything consistent.
        self.dV = reduce(lambda x, y: x * y, self.gridSize)

        self.regul = regul 
        self.cellStates = sim.cellStates


    def setBiophysics(self, biophysics):
        self.biophys = biophysics

    def setRegulator(self, regul):
        self.regul = regul
    
    def flattenIdx(self, idx):
        print "idx = " + str(idx)
        print "flat idx = " + str(idx[2] + idx[1]*self.gridDim[3] + idx[0]*self.gridDim[2]*self.gridDim[3])
        return idx[2] + idx[1]*self.gridDim[3] + idx[0]*self.gridDim[2]*self.gridDim[3]

    def idxFromPos(self, p):
        return ( floor((p[0] - self.gridOrig[0]) / self.gridSize[0]), \
            floor((p[1] - self.gridOrig[1]) / self.gridSize[1]), \
            floor((p[2] - self.gridOrig[2]) / self.gridSize[2]) )

    def trilinearWeights(self, p):
        x = (p[0] - self.gridOrig[0]) / self.gridSize[0]
        dx = x - floor(x)
        y = (p[1] - self.gridOrig[1]) / self.gridSize[1]
        dy = y - floor(y)
        z = (p[2] - self.gridOrig[2]) / self.gridSize[2]
        dz = z - floor(z)
        w = numpy.zeros((2,2,2))
        
        # Calculate weights, giving zero if grid point is off grid
        w[0,0,0] = (1-dx)*(1-dy)*(1-dz) if x>=0 and x<self.gridDim[1] and y>=0 and y<self.gridDim[2] and z>=0 and z<self.gridDim[3] else 0.0
        w[1,0,0] = dx*(1-dy)*(1-dz) if x>=-1 and x<self.gridDim[1]-1 and y>=0 and y<self.gridDim[2] and z>=0 and z<self.gridDim[3] else 0.0
        w[0,1,0] = (1-dx)*dy*(1-dz) if x>=0 and x<self.gridDim[1] and y>=-1 and y<self.gridDim[2]-1 and z>=0 and z<self.gridDim[3] else 0.0
        w[1,1,0] = dx*dy*(1-dz) if x>=-1 and x<self.gridDim[1]-1 and y>=-1 and y<self.gridDim[2]-1 and  z>=0 and z<self.gridDim[3] else 0.0
        w[0,0,1] = (1-dx)*(1-dy)*dz if x>=0 and x<self.gridDim[1] and y>=0 and y<self.gridDim[2] and z>=-1 and z<self.gridDim[3]-1 else 0.0
        w[1,0,1] = dx*(1-dy)*dz if x>=-1 and x<self.gridDim[1]-1 and y>=0 and y<self.gridDim[2] and z>=-1 and z<self.gridDim[3]-1 else 0.0
        w[0,1,1] = (1-dx)*dy*dz if x>=0 and x<self.gridDim[1] and y>=-1 and y<self.gridDim[2]-1 and z>=-1 and z<self.gridDim[3]-1 else 0.0
        w[1,1,1] = dx*dy*dz if x>=-1 and x<self.gridDim[1]-1 and y>=-1 and y<self.gridDim[2]-1 and z>=-1 and z<self.gridDim[3]-1 else 0.0
        print "w = "+str(w.reshape(8))
        return w 

    def dataLen(self):
        return self.gridDataLen 

    def addCell(self, cellState):
        pass

    def transportRates(self, signalRates, signalLevels, boundcond='constant', mode='normal'):
        # Compute diffusion term, laplacian of grid levels in signalLevels, 
        # write into signalRates
        #
        # mode='greens' - do not use initLevels as these don't apply!
        signalRatesView = signalRates.reshape(self.gridDim)
        signalLevelsView = signalLevels.reshape(self.gridDim)
        advKernel = numpy.zeros((3,3,3))
        advKernel[:,1,1] = [-0.5,0,0.5]
        for s in range(self.nSignals):
            if boundcond=='constant' and self.initLevels and mode!='greens':
                boundval = self.initLevels[s]
            else:
                boundval = 0.0 
            if self.advRates:
                # Adevction term = du/dx
                # Note: always use 'nearest' edge case, this gives central 
                # differences in middle, and forward/backward differences on edges
                convolve(signalLevelsView[s], advKernel*self.advRates[s], output=signalRatesView[s], mode='nearest')
                # Diffusion term = \del^2u
                # Use edge case from boundary conditions for diffusion
                signalRatesView[s] += laplace(signalLevelsView[s], None, mode=boundcond, cval=boundval) * self.diffRates[s]
            else:
                signalRatesView[s] = laplace(signalLevelsView[s], None, mode=boundcond, cval=boundval) * self.diffRates[s]

 
    def interpAddToGrid(self, pos, delta, grid):
        pidx = numpy.array(self.idxFromPos(pos))
        w = self.trilinearWeights(pos)
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    idx = tuple(pidx + numpy.array((i,j,k)))
                    for s in range(self.nSignals):
                        grid[(s,)+idx] += w[i,j,k]*delta[s]

    def cellProdRates(self, signalRates, cellState, cellRates):
        signalRatesView = signalRates.reshape(self.gridDim)
        self.interpAddToGrid(cellState.pos, cellRates, signalRatesView)
#        pidx = numpy.array(self.idxFromPos(cellState.pos))
#        w = self.trilinearWeights(cellState.pos)
#        for i in range(2):
#            for j in range(2):
#                for k in range(2):
#                    idx = tuple(pidx + numpy.array((i,j,k)))
#                    for s in range(self.nSignals):
#                        signalRatesView[(s,)+idx] += w[i,j,k]*cellRates[s]

    def signals(self, cellState, signalLevels):
        # Gets signal level at nearest neighbour grid point to cell position
        # TO DO: interpolation
        signalLevelsView = signalLevels.reshape(self.gridDim)
        pidx = self.idxFromPos(cellState.pos)
        print "flatidx = %i"%(self.flattenIdx(pidx))
        w = self.trilinearWeights(cellState.pos)
        sigs = numpy.zeros((self.nSignals))
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    idx = tuple(pidx + numpy.array((i,j,k),dtype=int))
                    for s in range(self.nSignals):
                        sigs[s] += w[i,j,k]*signalLevelsView[(s,)+idx]
        return sigs

    def initSignalLevels(self, levels):
        # levels already contains result of previous time step. Since the
        # signals are on a fixed grid, no changes needed.
        pass

    def step(self, dt):
        # Nothing specific to this model. Integrator does everything.
        pass
