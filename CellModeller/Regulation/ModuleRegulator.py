import copy
import os.path
import sys
import imp

class ModuleRegulator:
    def __init__(self, sim, biophys=None, signalling=None):
        self.modName = sim.moduleName
        self.modStr = sim.moduleStr
        self.sim = sim
        self.cellStates = sim.cellStates
        self.biophys = biophys
        self.signal = signalling
        # Simulator is responsible for loading the model as a python module
        # This class uses the module imported by Simulator
        self.module = sim.module 

    def addCell(self, cellState, **kwargs):
        self.module.init(cellState, **kwargs)

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
        csv = list(self.cellStates.values())
        nCells = len(csv)
        for i in range(nCells):
            levels[i,:] = csv[i].species

    def step(self, dt=0):
        try:
            self.module.update(self.cellStates)
        except Exception as e:
            print("Problem with regulation module " + self.modName)
            print(e)

    def divide(self, pState, d1State, d2State):
        # Call the module's optional divide function
        divfunc = getattr(self.module, "divide", None)
        if callable(divfunc):
            divfunc(pState, d1State, d2State)


