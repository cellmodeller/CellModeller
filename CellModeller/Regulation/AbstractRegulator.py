class AbstractRegulator:

    """Models the genetic regulation of cells.  This class exists only
    to document the required interface.  (There is no need to inherit
    from it.)
    """

    def __init__(self, simulator):
        """Reference to Simulator required.

        simulator -- The Simulator object containing CellStates.
        Should be used to access current cell states and other
        simulation parameters
         """
        raise NotImplementedError()

    def update(self, dt):
        """Step the simulation forward.

        dt -- Time to step forward.
        """
        raise NotImplementedError()

    def setBiophysics(self, biophysics):
        """Set the biophysical model to be used by this regulator.

        biophysics -- A biophysical model.
        """
        raise NotImplementedError()

    def setSignalling(self, signalling):
        """Set the signalling model to be used by this regulator.

        signalling -- A signalling model.
        """
        raise NotImplementedError()

    def addCell(self, cellState):
        """Add a cell to the model.

        cellState -- The state of the new cell.  Id should be unique.
        """
        raise NotImplementedError()

    def parameters(self, cellState):
        """Return the parameters for the biophysical model of a cell.

        cellState -- The state of the cell to return the parameters of.
        """
        raise NotImplementedError()

    def speciesRates(self, cellState, speciesLevels, signalLevels):
        """Return rates of species production (dydt) by a cell,
        given species and signal levels (y) supplied. I.e. dydt(y)

        Must return rates in a numpy array

        cellState -- The state of the cell to return rates of.
        speciesLevels -- Levels of species in the cell.
        signalLevels -- Levels of signals in the vicinity of the cell.
        """
        raise NotImplementedError()

    def signalRates(self, cellState, speciesLevels, signalLevels):
        """Return rates of signal production (dydt) by a cell,
        given species and signal levels (y) supplied. I.e. dydt(y)

        Must return rates in a numpy array

        cellState -- The state of the cell to return rates of.
        speciesLevels -- Levels of species in the cell.
        signalLevels -- Levels of signals in the vicinity of the cell.
        """
        raise NotImplementedError()

    def initSpeciesLevels(levels):
        """Set the initial species levels for time step to current levels.
        Write directly into levels array.

        levels -- the array stored by integrator
        """
        raise NotImplementedError()

    def isDividing(self, cellState):
        """Return True if this cell is dividing.

        cellState -- The state of the cell to return the division status of.
        """
        raise NotImplementedError()

    def divide(self, parent_id, daughter_id1, daughter_id2):
        """Divide a cell.

        parent_id -- The id of the cell to divide.
        daughter_id1 -- The id of the first daughter cell.  Should be unique.
        daughter_id2 -- The id of the second daughter cell.  Should be unique.
        """
        raise NotImplementedError()

