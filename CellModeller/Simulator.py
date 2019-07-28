from CellState import CellState
import copy
import pyopencl as cl
import sys
import os
import cPickle
import csv
import numpy
import inspect
import imp
import ConfigParser
import importlib

class Simulator:
    """
This class is in charge of running the simulation, creating the various models
and stepping them forward in time. It is the control interface for the gui
or script that is running the simulation.

Stores a map from cell_id to CellState, which stores the current simulation
state of each cell.

Constructed on a user-defined python file. This file implements a
function setup(Simulator, Gui) that constructs the requiredx modules
(Regulator, Signalling, Integrator), and calls Simulator.init(). It
can also create Renderers and add them by calling
Simulator.addRenderer(renderer) so that the simulation can be
visualised.
"""

    ## Construct an empty simulator object. This object will not be able to
    # do anything yet unti we use 'init' method to specify the models for
    # physical interaction, genetic circuit, diffusion and integrator.
    def __init__(self, \
                    moduleName, \
                    dt, \
                    pickleSteps=50, \
                    outputDirName=None, \
                    moduleStr=None, \
                    saveOutput=False, \
                    clPlatformNum=0, \
                    clDeviceNum=0, \
                    is_gui=False):
        # Is this simulator running in a gui?
        self.is_gui = is_gui

        # No models specified yet
        self.reg = None
        self.phys = None
        self.sig = None
        self.integ = None
        self.pickleSteps = pickleSteps

        # No cells yet, initialise indices and empty lists/dicts, zero counters
        self._next_id = 1
        self._next_idx = 0
        self.idToIdx = {}
        self.idxToId = {}
        self.cellStates = {}
        self.renderers = []
        self.stepNum = 0
        self.lineage = {}

        # Time step
        self.dt = dt

        if "CMPATH" in os.environ:
            self.cfg_file = os.path.join(os.environ["CMPATH"], 'CMconfig.cfg')
        else:
            self.cfg_file = 'CellModeller/CMconfig.cfg'
        if not self.init_cl(platnum=clPlatformNum, devnum=clDeviceNum):
            print "Couldn't initialise OpenCL context"
            return

        # Two ways to specify a module (model):
        self.moduleName = moduleName # Import via standard python
        self.moduleStr = moduleStr # Import stored python code string (from a pickle usually)
        if self.moduleStr:
            print "Importing model %s from string"%(self.moduleName)
            self.module = imp.new_module(moduleName)
            exec moduleStr in self.module.__dict__
        else:
            # In case the moduleName is a path to a python file:
            # Get path and file name
            (path,name) = os.path.split(self.moduleName)
            # Append path to PYTHONPATH, if no path do nothing
            if path:
                if path not in sys.path:
                    sys.path.append(path)
            # Remove .py extension if present
            self.moduleName = str(name).split('.')[0]
            print "Importing model %s"%(self.moduleName)
            if self.moduleName in sys.modules:
                self.module = sys.modules[self.moduleName]
                reload(self.module)
            else:
                self.module = __import__(self.moduleName, globals(), locals(), [], -1)
            

        # TJR: What is this invar thing? I have never seen this used...
        #setup the simulation here:
        #if invar:
        #    self.module.setup(self, invar)
        #else:

        # Set up the data output directory
        self.dataOutputInitialised=False
        self.outputDirName = outputDirName
        self.setSaveOutput(saveOutput)
        '''
        self.saveOutput = saveOutput
        if self.saveOutput:
            self.outputSteps = outputSteps
            self.init_data_output(outputFileDir)
        '''
        
        # Call the user-defined setup function on ourself
        self.module.setup(self)

    def setSaveOutput(self, save):
        self.saveOutput = save
        if save and (not self.dataOutputInitialised):
            self.init_data_output()

    def init_data_output(self):
        import time
        startTime = time.localtime()
        outputFileRoot = self.outputDirName if self.outputDirName else self.moduleName + '-' + time.strftime('%y-%m-%d-%H-%M', startTime)
        self.outputDirPath = os.path.join('data', outputFileRoot)
        if 'CMPATH' in os.environ:
            self.outputDirPath = os.path.join(os.environ["CMPATH"], self.outputDirPath)

        # Add a number to end of dir name if it already exists 
        label = 2
        while os.path.exists(self.outputDirPath):
            if label>2:
                self.outputDirPath = self.outputDirPath[:-2]+"_"+str(label)
            else:
                self.outputDirPath = self.outputDirPath+"_"+str(label)
            label+=1
        os.mkdir(self.outputDirPath)

        # write a copy of the model into the dir (for reference), 
        # this goes in the pickle too (and gets loaded when a pickle is loaded)
        if self.moduleStr:
            self.moduleOutput = self.moduleStr
        else:
            self.moduleOutput = inspect.getsource(self.module)
        open(os.path.join(self.outputDirPath, self.moduleName), 'w').write(self.moduleOutput)

        self.dataOutputInitialised=True


    ## Get an id for the next cell to be created
    def next_id(self):
        id = self._next_id
        self._next_id += 1
        return id

    ## Get the index (into flat arrays) of the next cell to be created
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

    ## Specify models to be used by simulator object. The four inputs are
    # 'phys' = physical model of cell iteractions
    # 'reg' = regulatory model of biochemical circuit in the cell
    # 'sig' = signaling model of intercellular chemical reaction diffusion.
    # 'integ' = integrator

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


    ## Set up the OpenCL contex, the configuration is set up the first time, and is saved in the config file
    def init_cl(self, platnum, devnum):
        # Check that specified platform exists
        platforms = cl.get_platforms()
        if len(platforms)<=platnum:
            print "Specified OpenCL platform number (%d) does not exist."
            print "Options are:"
            for p in range(len(platforms)):
                print "%d: %s"%(p, str(platforms[p])) 
            return False
        else:
            platform = platforms[platnum]

        # Check that specified device exists on that platform
        devices = platforms[platnum].get_devices()
        if len(devices)<=devnum:
            print "Specified OpenCL device number (%d) does not exist on platform %s."%(devnum,platform)
            print "Options are:"
            for d in range(len(devices)):
                print "%d: %s"%(d, str(devices[d])) 
            return False
        else:
            device = devices[devnum]

        # Create a context and queue
        self.CLContext = cl.Context(properties=[(cl.context_properties.PLATFORM, platform)],
                                          devices=[device])
        self.CLQueue = cl.CommandQueue(self.CLContext)
        print "Set up OpenCL context:"
        print "  Platform: %s"%(str(platform.name))
        print "  Device: %s"%(str(device.name))
        return True
        
    ## Get the OpenCL context and queue for running kernels 
    def getOpenCL(self):
        return (self.CLContext, self.CLQueue)

    ## set cell states from a given dict
    def setCellStates(self, cellStates):
        #Set cell states, e.g. from pickle file
        #this sets them on the card too
        self.cellStates = {}
        self.cellStates = cellStates
        idx_map = {}
        id_map = {}
        idmax = 0
        for id,state in cellStates.iteritems():
            idx_map[state.id] = state.idx
            id_map[state.idx] = state.id
            if id>idmax:
                idmax=id
        self.idToIdx = idx_map
        self.idxToId = id_map
        self._next_id = idmax+1
        self._next_idx = len(cellStates)
        self.reg.cellStates = cellStates
        self.phys.load_from_cellstates(cellStates)
    
    ## Add a graphics renderer - this should not be here --> GUI
    def addRenderer(self, renderer):
        self.renderers.append(renderer)
        
    ## Reset the simulation back to initial conditions
    def reset(self):
        # Delete existing models
        if self.phys:
            del self.phys
        if self.sig:
            del self.sig
        if self.integ:
            del self.integ
        if self.reg:
            del self.reg

        if not self.moduleStr: 
            #This will take up any changes made in the model file
            reload(self.module)
        else:
            # TJR: Module loaded from pickle, cannot reset?
            pass

        # Lose old cell states
        self.cellStates = {}
        # Recreate models via module setup
        self.module.setup(self)


    # Divide a cell to two daughter cells
    def divide(self, pState):
        pState.divideFlag = False
        pid = pState.id
        d1id = self.next_id()
        d2id = self.next_id()
        d1State = copy.deepcopy(pState)
        d2State = copy.deepcopy(pState)
        d1State.id = d1id
        d2State.id = d2id

        #reset cell ages
        d1State.cellAge = 0
        d2State.cellAge = 0
        #inherit effGrowth
        d1State.effGrowth = pState.effGrowth
        d2State.effGrowth = pState.effGrowth
        
        self.lineage[d1id] = pid
        self.lineage[d2id] = pid

        # Update CellState map
        self.cellStates[d1id] = d1State
        self.cellStates[d2id] = d2State
        del self.cellStates[pid]

        # Update indexing, reuse parent index for d1
        d1State.idx = pState.idx
        self.idToIdx[d1id] = pState.idx
        self.idxToId[pState.idx] = d1id
        d2State.idx = self.next_idx()
        self.idToIdx[d2id] = d2State.idx
        self.idxToId[d2State.idx] = d2id
        del self.idToIdx[pid]

        # Divide the cell in each model
        asymm = getattr(pState, 'asymm', [1, 1])
        self.phys.divide(pState, d1State, d2State, f1=asymm[0], f2=asymm[1])
        if self.integ:
            self.integ.divide(pState, d1State, d2State)
        self.reg.divide(pState, d1State, d2State)

    ## Add a new cell to the simulator
    def addCell(self, cellType=0, cellAdh=0, length=3.5, **kwargs):
        cid = self.next_id()
        cs = CellState(cid)
        cs.length = length
        cs.cellType = cellType
        cs.cellAdh = cellAdh
        cs.idx = self.next_idx()
        self.idToIdx[cid] = cs.idx
        self.idxToId[cs.idx] = cid
        self.cellStates[cid] = cs
        if self.integ:
            self.integ.addCell(cs)
        self.reg.addCell(cs)
        if self.sig:
            self.sig.addCell(cs)
        self.phys.addCell(cs, **kwargs)

    #---
    # Some functions to modify existing cells (e.g. from GUI)
    # Eventually prob better to have a generic editCell() that deals with this stuff
    #
    def moveCell(self, cid, delta_pos):
        if self.cellStates.has_key(cid):
            self.phys.moveCell(self.cellStates[cid], delta_pos)

    ## Proceed to the next simulation step
    # This method is where objects phys, reg, sig and integ are called
    def step(self):
        self.reg.step(self.dt)
        states = dict(self.cellStates)
        for (cid,state) in states.items():
            if state.divideFlag:
                self.divide(state) #neighbours no longer current

        self.phys.set_cells()
        while not self.phys.step(self.dt): #neighbours are current here
            pass
        if self.sig:
            self.sig.step(self.dt)
        if self.integ:
            self.integ.step(self.dt)

        if self.saveOutput and self.stepNum%self.pickleSteps==0:
            self.writePickle()

        self.stepNum += 1
        return True


    ## Import cells to the simulator from csv file. The file contains a list of 7-coordinates {pos,dir,len} (comma delimited) of each cell - also, there should be no cells around - ie run this from an empty model instead of addcell
    def importCells_file(self, filename):
        f=open(filename, 'rU')
        list=csv.reader(f,delimiter=',')
        for row in list:
            cpos = [float(row[0]),float(row[1]),float(row[2])]
            cdir = [float(row[3]),float(row[4]),float(row[5])]
            clen = float(row[6]) #radius should be removed from this in the analysis
            ndir = cdir/numpy.linalg.norm(cdir) #normalize cell dir just in case
            #this should probably also check for overlaps
            self.addCell(pos=tuple(cpos), dir=tuple(ndir), length=clen)

    ## Write current simulation state to an output file
    def writePickle(self, csv=False):
        filename = os.path.join(self.outputDirPath, 'step-%05i.pickle' % self.stepNum)
        outfile = open(filename, 'wb')
        data = {}
        data['cellStates'] = self.cellStates
        data['stepNum'] = self.stepNum
        data['lineage'] = self.lineage
        data['moduleStr'] = self.moduleOutput
        data['moduleName'] = self.moduleName
        if self.integ:
            print("Writing new pickle format")
            data['specData'] = self.integ.levels
            data['sigGrid'] = self.integ.signalLevel
            data['sigGridOrig'] = self.sig.gridOrig
            data['sigGridDim'] = self.sig.gridDim
            data['sigGridSize'] = self.sig.gridSize
        if self.sig:
            data['sigData'] = self.integ.cellSigLevels
            data['sigGrid'] = self.integ.signalLevel
        cPickle.dump(data, outfile, protocol=-1)
        #output csv file with cell pos,dir,len - sig?

    # Populate simulation from saved data pickle
    def loadFromPickle(self, data):
        self.setCellStates(data['cellStates'])
        self.lineage = data['lineage']
        self.stepNum = data['stepNum']
        idx_map = {}
        id_map = {}
        idmax = 0
        for id,state in data['cellStates'].iteritems():
            idx_map[state.id] = state.idx
            id_map[state.idx] = state.id
            if id>idmax:
                idmax=id
        self.idToIdx = idx_map
        self.idxToId = id_map
        self._next_id = idmax+1
        self._next_idx = len(data['cellStates'])
        if self.integ:
            if data.has_key('sigData'):
                self.integ.setLevels(data['specData'],data['sigData'])
            elif data.has_key('specData'):
                self.integ.setLevels(data['specData'])


