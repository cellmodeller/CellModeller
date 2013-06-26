import numpy
import scipy.integrate.odepack

class ScipyODEIntegrator:
    def __init__(self, sim, nSignals, nSpecies, maxCells, sig=None, fixedSignalLen=False, regul=None):
        self.regul = regul
        self.cellStates = sim.cellStates
        self.nCells = len(self.cellStates)

        self.nSpecies = nSpecies
        self.nSignals = nSignals
        self.maxCells = maxCells
        self.fixedSignalLen = fixedSignalLen

        # scipy ode solver
        self.ode = scipy.integrate.ode(self.dydt)\
            .set_integrator('vode', method='bdf', with_jacobian=False, nsteps=10000, rtol=1e-3)
         
        # The signalling model
        self.signalling = sig

             # Data organisation
        if self.signalling:
            self.signalDataLen = self.signalling.dataLen()
        else:
            self.signalDataLen = 0
        self.maxSpecDataLen = self.maxCells*nSpecies
        if fixedSignalLen:
            # no need to scale up signal storage
            storageLen = self.maxSpecDataLen + self.signalDataLen
        else:
            # allow for storage of maxCells
            storageLen = self.maxSpecDataLen + self.nSignals*self.maxCells

        # These arrays store the level and rate of signals and species
        # in a contiguous form. The first part is the signals,
        # then the cell species
        # To avoid reallocation, create enough space for maxCells

        self.levels = numpy.zeros(storageLen)
        self.rates = numpy.zeros(storageLen)
        self.makeViews()

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
    
    def dydt(self, t, y):
        # Make some views on the current levels, y
        signalLvlView = y[0:self.signalDataLen]
        speciesLvlView = y[self.signalDataLen:self.dataLen].reshape(self.nCells,self.nSpecies)

        # The following writes the rates into the correct parts of self.rates,
        # using the views created in constructor

        # Get transport term from signalling model
        if self.signalling:
            self.signalling.transportRates(self.signalRate, signalLvlView)

        # Loop over cells to get rates
        states = self.cellStates
        for (id,c) in states.items():
            if self.signalling:
                cellSignals = self.signalling.signals(c, signalLvlView)
            else:
                cellSignals = []
            idx = c.idx
            self.specRate[idx,:] = self.regul.speciesRates(c, speciesLvlView[idx,:], cellSignals)
            cellRates = self.regul.signalRates(c, speciesLvlView[idx,:], cellSignals)
            if self.signalling:
                self.signalling.cellProdRates(self.signalRate, c, cellRates)

        # Return the rate array (up to the actual number of cells)
        return self.rates[0:self.dataLen]

    def step(self, dt):
        #self.nCells = len(self.cellStates)
        # Check we have enough space allocated
        try:
            s = self.specLevel[self.nCells-1]
        except IndexError:
            # Could resize here, then would have to rebuild views
            print "Number of cells exceeded " \
                + self.__class__.__name__ \
                + "::maxCells (" + self.maxCells + ")"
           
        if self.signalling and not self.fixedSignalLen:
            # Check if size of signalling part changed (new cells)
            signalDataLen = self.signalling.dataLen()
            if signalDataLen>self.signalDataLen:
                # Change the views into signal/species parts of data
                self.signalDataLen = signalDataLen
                self.makeViews()

        self.dataLen = self.signalDataLen + self.nCells*self.nSpecies
        # Set initial conditions for this step
        #self.regul.initSpeciesLevels(self.specLevel[0:self.nCells,:])
        #if self.signalling:
        #    self.signalling.initSignalLevels(self.signalLevel)
        self.ode.set_initial_value(self.levels[0:self.dataLen], 0)
        # Integrate the ode
        self.ode.integrate(dt)
        self.levels[0:self.dataLen] = self.ode.y

        # Put the final signal levels into the cell states
        states = self.cellStates
        for (id,c) in states.items():
            if self.signalling:
                c.signals = self.signalling.signals(c, self.signalLevel)

            
 
