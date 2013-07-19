class AbstractIntegrator:

    """Integrates signalling and genetic regulation. Takes rates from Regulator
    and combines with signal transport, stepping the system forward in time.

     This class exists only to document the required interface.  (There is no need to inherit
    from it.)
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
        """Step the model forward by dt, updating the CellState contents.

               dt -- Time to step forward.
        """
        raise NotImplementedError()

    def addCell(self, cellState):
        """ Add cell to the model

        cellState -- state of new cell. Id should be unique.
        """

    def setRegulator(self, regulator):
        """Set the regulator object used by this integrator 

        regulator -- A regulator.
        """
        raise NotImplementedError()

    def setSignalling(self, signalling):
        """Set the signalling model object used by this integrator 

        signalling -- A signalling.
        """
        raise NotImplementedError()

