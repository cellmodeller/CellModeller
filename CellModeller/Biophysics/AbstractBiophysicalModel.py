class AbstractBiophysicalModel:

    """Biophysical model of cells.  This class exists only to document
    the required interface.  (There is no need to inherit from it.)
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

    def step(self, dt):
        """Step the simulation forward. Also update geometry information
            in CellState for each cell.

        dt -- Time to step forward.
        """
        raise NotImplementedError()

    def setRegulator(self, regulator):
        """Specify the Regulator object to be used by the model.

        regulator -- A regulator.
        """
        raise NotImplementedError()

    def addCell(self, cellState, *args, **keywords):
        """Add a cell to the simulation.

        Implementations may require additional keyword arguments to
        specify e.g. cell location or size.

        cellState -- The state of the added cell.  Id should be unique.
        """
        raise NotImplementedError()

    def hasNeighbours(self):
        """ Does this model implement neighbours? Return True/False
        """

    def divide(self, parentState, daughter1State, daughter2State, *args, **kw):
        """Divide a cell.

        parentState -- The state of the cell to divide.
        daughter1State -- The state of the first daughter cell.  Should be unique.
        daughter2State -- The state of the second daughter cell.  Should be unique.

        Additional keywords args for model specific behaviour.
        """
        raise NotImplementedError()

