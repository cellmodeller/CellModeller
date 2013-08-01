from CellState import CellState
import copy
import pyopencl as cl
import sys
import os
import cPickle

class Simulator:
    """
This class is in charge of running the simulation, creating the various models
and stepping them forward in time. It is the control interface for the gui
or script that is running the simulation.

Stores a map from cell_id to CellState, which stores the current simulation
state of each cell.

Constructed on a user-defined python file. This file implements a
function setup(Simulator, Gui) that constructs the required modules
(Regulator, Signalling, Integrator), and calls Simulator.init(). It
can also create Renderers and add them by calling
Simulator.addRenderer(renderer) so that the simulation can be
visualised.
"""


    def __init__(self, moduleName, dt, pickleSteps=50, pickleFileRoot=False):
        self.dt = dt
        self._next_id = 1
        self._next_idx = 0
        self.idToIdx = {}
        # Map from id to cell state
        self.cellStates = {}
        self.reg = None
        self.phys = None
        self.sig = None
        self.integ = None
        self.renderers = []
        self.stepNum = 0
        self.savePickle =True 
        self.pickleSteps = pickleSteps

        self.lineage = {}

        self.init_cl()

        if moduleName:
            self.moduleName = moduleName
            self.module = __import__(self.moduleName, globals(), locals(), [], -1)
            self.module.setup(self)
            import time
            self.startTime = time.localtime()
            self.pickleFileRoot = pickleFileRoot if pickleFileRoot else self.moduleName + '-' + time.strftime('%H-%M-%d-%m-%y', self.startTime)
            self.pickleDir = os.path.join('data', self.pickleFileRoot)
            os.mkdir(self.pickleDir) # raises OSError if dir already exists
            # write a copy of the model into the dir (for reference)
            self.moduleStr = open(self.module.__file__, 'rU').read()
            open(os.path.join(self.pickleDir, self.moduleName), 'w').write(self.moduleStr)


    def next_id(self):
        id = self._next_id
        self._next_id += 1
        return id

    def next_idx(self):
        idx = self._next_idx
        self._next_idx += 1
        return idx

    # Currently, the user-defined regulation module creates the
    # biophysics, regulation, and signalling objects in a function
    # setup().
    #
    # We pass in the empty simulator object, ie. setup(sim)
    # and have the user-defined func initialise the 3 modules
    def init(self, phys, reg, sig, integ):
        self.phys = phys
        self.reg = reg
        self.sig = sig

        self.phys.setRegulator(reg)

        self.reg.setBiophysics(phys)

        if integ:
            self.integ = integ
            self.integ.setRegulator(reg)

        if self.sig:
            self.sig.setBiophysics(phys)
            self.sig.setRegulator(reg)
            self.integ.setSignalling(sig)
            self.reg.setSignalling(sig)


    def init_cl(self):
        """Set up the OpenCL context."""
        platform = cl.get_platforms()[0]
        if sys.platform == 'darwin':
            self.CLContext = cl.Context(devices=[platform.get_devices()[0]])
        else:
            #try:
            #    self.CLContext 
            #    = cl.Context(properties=[(cl.context_properties.PLATFORM, 
            #    platform)])
            #except:
            self.CLContext = cl.Context(properties=[(cl.context_properties.PLATFORM, platform)],
                                          devices=[platform.get_devices()[0]])
        self.CLQueue = cl.CommandQueue(self.CLContext)
        print platform.get_devices()[0]
        print (platform.get_devices()[0]).get_info(cl.device_info.DRIVER_VERSION)

    def getOpenCL(self):
        return (self.CLContext, self.CLQueue)

    def setCellStates(self, cellStates):
        # Set cell states, e.g. from pickle file 
        self.cellStates = cellStates

    def addRenderer(self, renderer):
        self.renderers.append(renderer)

    def reset(self):
        # Delete the models, they might be holding up memory/GPU resources?
        if self.phys:
            del self.phys
        if self.sig:
            del self.sig
        if self.integ:
            del self.integ
        if self.reg:
            del self.reg
        # Lose old cell states
        self.cellStates = {}
        # Recreate models via module setup
        self.module.setup(self)
        #self.phys.reset()
        #if self.sig:
        #    self.sig.reset()
        #self.integ.reset()

    def divide(self, pState):
        pState.divideFlag = False
        pid = pState.id
        d1id = self.next_id()
        d2id = self.next_id()
        d1State = copy.deepcopy(pState)
        d2State = copy.deepcopy(pState)
        d1State.id = d1id
        d2State.id = d2id

        self.lineage[d1id] = pid
        self.lineage[d2id] = pid

        # Update CellState map
        self.cellStates[d1id] = d1State
        self.cellStates[d2id] = d2State
        del self.cellStates[pid]

        # Update indexing, reuse parent index for d1
        d1State.idx = pState.idx
        self.idToIdx[d1id] = pState.idx
        d2State.idx = self.next_idx()
        self.idToIdx[d2id] = d2State.idx
        del self.idToIdx[pid]

        # Divide the cell in each model
        asymm = getattr(pState, 'asymm', [1, 1])
        self.phys.divide(pState, d1State, d2State, f1=asymm[0], f2=asymm[1])
        if self.integ:
            self.integ.divide(pState, d1State, d2State)
        self.reg.divide(pState, d1State, d2State)


    def addCell(self, cellType=0, **kwargs):
        cid = self.next_id()
        cs = CellState(cid)
        cs.idx = self.next_idx()
        cs.cellType = cellType
        self.idToIdx[cid] = cs.idx
        self.cellStates[cid] = cs
        if self.integ:
            self.integ.addCell(cs)
        self.reg.addCell(cs)
        if self.sig:
            self.sig.addCell(cs)
        self.phys.addCell(cs, **kwargs)


    def step(self):
        self.reg.step(self.dt)
        if self.sig:
            self.sig.step(self.dt)
        self.phys.step(self.dt)
        if self.integ:
            self.integ.step(self.dt)

        states = dict(self.cellStates)
        for (cid,state) in states.items():
            if state.divideFlag:
                self.divide(state)

        if self.savePickle and self.stepNum%self.pickleSteps==0:
            self.writePickle()

        self.stepNum += 1


    def writePickle(self):
        filename = os.path.join(self.pickleDir, 'step-%05i.pickle' % self.stepNum)
        outfile = open(filename, 'wb')
        if self.integ and self.sig:
            sigData = (self.sig.gridSize, self.sig.gridOrig, self.integ.gridDim, self.integ.signalLevel.reshape(self.integ.gridDim))
            data = (self.cellStates, sigData, self.lineage)
        else:
            data = (self.cellStates, self.lineage)
        cPickle.dump(data, outfile, protocol=-1)
