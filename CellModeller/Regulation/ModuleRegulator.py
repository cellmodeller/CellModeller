import copy
import os.path
import sys

class ModuleRegulator:
    def __init__(self, sim, modName, biophys=None, signalling=None):
        self.modName = modName
        self.sim = sim
        self.cellStates = sim.cellStates
        self.biophys = biophys
        self.signal = signalling
        self.module = None
        self.reset()

    def addCell(self, cellState):
        self.module.init(cellState, sim=self.sim)

    def importModule(self):
        # modName may be a full path to file
        # if so, add the path, and split out
        # the module name
        (path,name) = os.path.split(self.modName)
        mod = str(name).split('.')[0]
        if path:
            if path not in sys.path:
                sys.path.append(path)
        self.module = __import__(mod, globals(), locals(), [], -1)

    def reset(self, modName=None):
        print(modName)
        if (modName and modName!=self.modName) or not self.module:
            if modName:
                self.modName = modName
            self.importModule()
        else:
            reload(self.module)

        self.nSpecies = self.module.numSpecies()
        self.nSignals = self.module.numSignals()

    def setSignalling(self, signal):
        self.signal = signal

    def setIntegrator(self, integ):
        self.integ = integ

    def setBiophysics(self, biophys):
        self.biophys = biophys

    def sigRateCL(self):
        return self.module.sigRateCL()

    def specRateCL(self):
        return self.module.specRateCL()

    def signalRates(self, cstate, speciesLevels, signalLevels):
        return self.module.signalRates(cstate, speciesLevels, signalLevels)

    def speciesRates(self, cstate, speciesLevels, signalLevels):
        return self.module.speciesRates(cstate, speciesLevels, signalLevels)

    def initSpeciesLevels(self, levels):
        csv = self.cellStates.values()
        nCells = len(csv)
        for i in range(nCells):
            levels[i,:] = csv[i].species

    def step_cell(self, dt, cellState):
        self.module.update_cell(dt, cellState, self.sim)

    def step(self, dt=0):
        #try:
        self.module.update(self.cellStates, self.sim)
        #except Exception as e:
        #    print "Problem with regulation module " + self.modName
        #    print e

    def divide(self, pState, d1State, d2State):
        # Call the module's optional divide function
        divfunc = getattr(self.module, "divide", None)
        if callable(divfunc):
            divfunc(pState, d1State, d2State)


