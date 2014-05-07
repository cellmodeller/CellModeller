import numpy

class NeighbourDiffusion:

    """Diffusion between neighbouring cells.
    """

    biophys = None
    regul = None

    def __init__(self, simulator, nSignals, D):
        self.sim = simulator
        self.nSignals = nSignals
        self.diffRates = numpy.array(D)
        self.nCells = 0

    def reset(self):
        """ Reset the model to its original condition, before any calls to step were made.
        """
        pass

    def dataLen(self):
        """ Returns the length of the data vector (array) that signals
        are stored/computed on
        """
        return self.nCells * self.nSignals

    def addCell(self, cellState):
        """ Add cell to the model

        cellState -- state of new cell. Id should be unique.
        """
        self.nCells +=1

    def setBiophysics(self, biophysics):
        """Set the biophysical model used by this signalling model.

        biophysics -- A biophysical model.
        """
        self.biophys = biophysics

    def setRegulator(self, regulator):
        """Set the regulator object used by this signalling model.

        regulator -- A regulator.
        """
        self.regul = regulator

    def signals(self, cellState, signalLevels):
        """Return the levels of signals near a cell, given the current
        levels in signalLevels array (1d array shape)

        cellState -- The state of the cell.
        """
        sigLvl = signalLevels.reshape(self.nCells,self.nSignals)
        idx = self.idToIdx[cellState.id]
        return sigLvl[idx,:]
        
    def initSignalLevels(self, levels):
        """Set the initial signal levels for time step to current levels.
        Write directly into levels array. The array will contain the 
        result of the previous time-step, so if there are no changes
        this function can be empty.

        levels -- the array stored by integrator
        """
        self.idToIdx = {}
        lvlView = levels.reshape(self.nCells, self.nSignals)
        for i in range(self.nCells):
            cell = self.sim.cellStates.values()[i]
            lvlView[i,:] = cell.signals
            self.idToIdx[cell.id] = i
        
    def transportRates(self, signalRates, signalLevels):
        """Write the rate of change of signals due to transport into signalRates, given the
        current signal levels in signalLevels. The two args are views into
        the data stored by the Integrator, and will be flat (1d). Reshape
        if necessary.

        signalRates -- view into numpy.array (1d) for writing result
        speciesRates -- view into numpy.array (1d) containing current signal levels
        """
        sigRates = signalRates.reshape(self.nCells, self.nSignals)
        sigLevels = signalLevels.reshape(self.nCells, self.nSignals)
        states = self.sim.cellStates.values()
        for i in range(self.nCells):
            cell = states[i]
            sigRates[i,:] = 0.0
            for n in cell.nbs:
                nbr = self.sim.cellStates[n]
                nidx = self.idToIdx[nbr.id]
                # This should be scaled by wall/contact area
                sigRates[i,:] += sigLevels[nidx,:]*self.diffRates
                sigRates[i,:] -= sigLevels[i,:]*self.diffRates

    def cellProdRates(self, signalRates, cellState, cellRates):
        """Add rate of change of signal due to cell production into signalRates

        signalRates -- array to write into
        cellState -- the cell
        cellRates -- the signal production rates
        """
        sigRates = signalRates.reshape(self.nCells,self.nSignals)
        idx = self.idToIdx[cellState.id]
        sigRates[idx,:] += cellRates

    def step(self, dt):
        self.nCells = len(self.sim.cellStates)
