class AbstractSignallingModel:

    """Signalling model.  This class exists only to document the
    required interface.  (There is no need to inherit from it.)
    """

   def __init__(self, simulator):
        """Reference to Simulator required.

        simulator -- The Simulator object containing CellStates.
        Should be used to access current cell states and other
        simulation parameters
        """
        raise NotImplementedError()

    def reset(self):
        """ Reset the model to its original condition, before any calls to step were made.
        """

    def dataLen(self):
        """ Returns the length of the data vector (array) that signals
        are stored/computed on
        """

    def addCell(self, cellState):
        """ Add cell to the model

        cellState -- state of new cell. Id should be unique.
        """

    def setBiophysics(self, biophysics):
        """Set the biophysical model used by this signalling model.

        biophysics -- A biophysical model.
        """
        raise NotImplementedError()

    def setRegulator(self, regulator):
        """Set the regulator object used by this signalling model.

        regulator -- A regulator.
        """
        raise NotImplementedError()

    def signals(self, cellState, signalLevels):
        """Return the levels of signals near a cell, given the current
        levels in signalLevels array (1d array shape)

        cellState -- The state of the cell.
        """
        raise NotImplementedError()

    def initSignalLevels(levels):
        """Set the initial signal levels for time step to current levels.
        Write directly into levels array. The array will contain the 
        result of the previous time-step, so if there are no changes
        this function can be empty.

        levels -- the array stored by integrator
        """
        raise NotImplementedError()

    def transportRates(self, signalRates, signalLevels):
        """Write the rate of change of signals due to transport into signalRates, given the
        current signal levels in signalLevels. The two args are views into
        the data stored by the Integrator, and will be flat (1d). Reshape
        if necessary.

        signalRates -- view into numpy.array (1d) for writing result
        speciesRates -- view into numpy.array (1d) containing current signal levels
        """
    def cellProdRates(self, signalRates, cellState, cellRates):
        """Add rate of change of signal due to cell production into signalRates

        signalRates -- array to write into
        cellState -- the cell
        cellRates -- the signal production rates
        """
        
    def step(self, dt):
        """ Do specific things for this model, for each time step
        """
