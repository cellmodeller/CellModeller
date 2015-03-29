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

    ## Construct an empty simulator object. This object will not be able to
    # do anything yet unti we use 'init' method to specify the models for
    # physical interaction, genetic circuit, diffusion and integrator.
    def __init__(self, moduleName, dt, invar=False, pickleSteps=50, pickleFileRoot=False, fromPickle=False):
        self.dt = dt
        self._next_id = 1
        self._next_idx = 0
        self.idToIdx = {}
        self.idxToId = {}
        # Map from id to cell state
        self.cellStates = {}
        self.reg = None
        self.phys = None
        self.sig = None
        self.integ = None
        self.renderers = []
        self.stepNum = 0
        self.savePickle = True
        self.pickleSteps = pickleSteps
        self.lineage = {}
        if fromPickle == False:
            self.fromPickle = False
        else:
            self.fromPickle = True

        if "CMPATH" in os.environ:
            self.cfg_file = os.path.join(os.environ["CMPATH"], 'CMconfig.cfg')
        else:
            self.cfg_file = 'CellModeller/CMconfig.cfg'
        self.init_cl()

        self.moduleName = moduleName
        print moduleName
        if fromPickle:
            self.moduleStr = fromPickle
            self.module = imp.new_module(moduleName)
            exec fromPickle in self.module.__dict__
        else:
            self.module = __import__(self.moduleName, globals(), locals(), [], -1)
        #setup the simulation here:
        if invar:
            self.module.setup(self, invar)
        else:
            self.module.setup(self)

        import time
        self.startTime = time.localtime()
        self.pickleFileRoot = pickleFileRoot if pickleFileRoot else self.moduleName + '-' + time.strftime('%y-%m-%d-%H-%M', self.startTime)
        self.pickleDir = os.path.join('data', self.pickleFileRoot)
        if "CMPATH" in os.environ:
            self.pickleDir = os.path.join(os.environ["CMPATH"], self.pickleDir)
        label = 2
        while os.path.exists(self.pickleDir):
            if label>2:
                self.pickleDir = self.pickleDir[:-2]+"_"+str(label)
            else:
                self.pickleDir = self.pickleDir+"_"+str(label)
            label+=1
        os.mkdir(self.pickleDir)
        # write a copy of the model into the dir (for reference), 
        # this goes in the pickle too (and gets loaded when a pickle is loaded)
        if not self.fromPickle:
            self.moduleStr = inspect.getsource(self.module)
        open(os.path.join(self.pickleDir, self.moduleName), 'w').write(self.moduleStr)


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
    # 'phy' = physical model of cell iteractions
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
    def init_cl(self):
        #if config file exists read and set everything
        config = ConfigParser.RawConfigParser()
        platform = cl.get_platforms()
        if os.path.isfile(self.cfg_file):
            config.read(self.cfg_file)
            platnum = int(config.get('Platform','platnum'))
            devnum = int(config.get('Device','devnum'))
        else:
            #select platform and device here and write to config file
            if len(platform) > 1:
                print "Select platform from the list:"
                for i in range(len(platform)):
                    print 'press '+str(i)+' for '+str(platform[i])
                platnum = int(input('Platform Number: '))
            else:
                platnum = 0
            if len(platform[platnum].get_devices())==0:
                print "No compatible device, check if your hardware is OpenCL compatible, and that OpenCL is correctly configured..." #this will now break
            elif len(platform[platnum].get_devices())==1:
                devnum = 0
            else:
                print "Select device from list"
                for i in range(len(platform[platnum].get_devices())):
                     print 'press '+str(i)+' for '+str(platform[platnum].get_devices()[i])
                devnum= int(input('Device number: '))
            #write settings to config file
            config.add_section('Platform')
            config.add_section('Device')
            config.set('Platform','platnum',platnum)
            config.set('Platform', 'Platform Name', platform[platnum])
            config.set('Device','devnum',devnum)
            config.set('Device', 'Device Name', platform[platnum].get_devices()[devnum])
            config.write(open(self.cfg_file, 'wb'))
        self.CLContext = cl.Context(properties=[(cl.context_properties.PLATFORM, platform[platnum])],
                                          devices=[platform[platnum].get_devices()[devnum]])
        self.CLQueue = cl.CommandQueue(self.CLContext)
        print (platform[platnum].get_devices()[devnum])
#.get_info(cl.device_info.DRIVER_VERSION)

    ## ??
    def getOpenCL(self):
        return (self.CLContext, self.CLQueue)

    ## set cell state

    def setCellStates(self, cellStates):
        #Set cell states, e.g. from pickle file
        #this sets them on the card too
        self.cellStates = {}
        self.cellStates = cellStates
        self.reg.cellStates = cellStates
        self.phys.load_from_cellstates(cellStates)
    
    ## ??
    def addRenderer(self, renderer):
        self.renderers.append(renderer)
        
    ## Delete the models, they might be holding up memory/GPU resources?
    def reset(self):
        if self.phys:
            del self.phys
        if self.sig:
            del self.sig
        if self.integ:
            del self.integ
        if self.reg:
            del self.reg
        if self.fromPickle == False: #This will take up any changes made in the model file
            reload(self.module)
        # Lose old cell states
        self.cellStates = {}
        # Recreate models via module setup
        self.module.setup(self)
        #self.phys.reset()
        #if self.sig:
        #    self.sig.reset()
        #self.integ.reset()


    # Divide a cell to two daughter cells
    def divide(self, pState):
        pState.divideFlag = False
        pid = pState.id
        d1id = self.next_id()
        d2id = self.next_id()
        d1State = copy.copy(pState)
        d2State = copy.copy(pState)
        d1State.id = d1id
        d2State.id = d2id

        #reset cell ages
        d1State.cellAge = 0
        d2State.cellAge = 0
        
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
    def addCell(self, cellType=0, cellAdh=0.0, length=3.5, **kwargs):
        cid = self.next_id()
        cs = CellState(cid)
        cs.length = length
        cs.oldLen = length
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
        return cid, cs

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


    ## Import cells from a file to the simulator from csv file. The file contains a list of 7-coordinates {pos,dir,len} (comma delimited) of each cell - also, there should be no cells around - ie run this from an empty model instead of addcell
    def importCells_file(self, filename):
        f=open(filename, 'rU')
        list=csv.reader(f,delimiter=',')
        for row in list:
            cpos = [float(row[0]),float(row[1]),float(row[2])]
            cdir = [float(row[3]),float(row[4]),float(row[5])]
            clen = float(row[6]) #radius should be removed from this in the analysis
            ndir = cdir/numpy.linalg.norm(cdir) #normalize cell dir just in case
            #this should probably also check for overlaps
            self.addCell(pos=tuple(cpos), dir=tuple(ndir), len=clen)

    ## Write current simulation state to an output file
    def writePickle(self, csv=False):
        filename = os.path.join(self.pickleDir, 'step-%05i.pickle' % self.stepNum)
        outfile = open(filename, 'wb')
        data = {}
        data['cellStates'] = self.cellStates
        data['stepNum'] = self.stepNum
        data['lineage'] = self.lineage
        data['moduleStr'] = self.moduleStr
        data['moduleName'] = self.moduleName
        if self.integ:
            data['specData'] = self.integ.levels
        if self.sig:
            data['sigData'] = self.integ.cellSigLevels
        cPickle.dump(data, outfile, protocol=-1)
        #output csv file with cell pos,dir,len - sig?

    def loadFromPickle(self, data):
        self.setCellStates(data['cellStates'])
        self.lineage = data['lineage']
        self.stepNum = data['stepNum']
        idx_map = {}
        idmax = 0
        for id,state in data['cellStates'].iteritems():
            idx_map[state.id] = state.idx
            if id>idmax:
                idmax=id
        self.idToIdx = idx_map
        self._next_id = idmax+1
        self._next_idx = len(data['cellStates'])
        if data.has_key('sigData'):
            self.integ.setLevels(data['specData'],data['sigData'])
        elif data.has_key('specData'):
            self.integ.setLevels(data['specData'])


