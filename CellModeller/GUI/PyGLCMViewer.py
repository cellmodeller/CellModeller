import PyQt4
from PyQt4 import QtCore, QtGui
from PyQt4.Qt import Qt
from PyQt4.QtCore import QObject, QTimer, pyqtSignal, pyqtSlot
from PyGLWidget import PyGLWidget
from OpenGL.GL import *
from OpenGL.GLU import *

from CellModeller.Regulation import ModuleRegulator
from CellModeller.Simulator import Simulator
from CellModeller.CellState import CellState
import os
import sys


class PyGLCMViewer(PyGLWidget):

    selectedCell = pyqtSignal(str)#CellState, name='selectedCell')
    selectedName = -1
    dt = 0.25 

    def __init__(self, parent = None):
        PyGLWidget.__init__(self,parent)
        self.animTimer = QTimer()
        self.animTimer.timeout.connect(self.animate)
        self.renderInfo = None
        self.sim= None
        self.modfile = None
        self.record = False
        self.set_radius(32)
        self.frameNo = 0

    def help(self):
        pass

    def setSimulator(self, sim):
        self.sim = sim

    @pyqtSlot(bool)
    def toggleRun(self, run):
        if run:
            self.animTimer.start(0)
        else:
            self.animTimer.stop()

    @pyqtSlot(bool)
    def toggleRecord(self, rec):
        self.record = rec
        self.sim.savePickle = rec

    @pyqtSlot()
    def reset(self):
        self.sim = Simulator(self.modname, self.dt)
        #if self.sim:
        #    self.sim.reset()
        self.frameNo = 0

    @pyqtSlot()
    def load(self):
        qs = QtGui.QFileDialog.getOpenFileName(self, 'Load Python module', '', '*.py')
        self.modfile = str(qs)
        self.loadFile(self.modfile)

    def loadFile(self, modstr):
        (path,name) = os.path.split(modstr)
        modname = str(name).split('.')[0]
        self.modname = modname
        sys.path.append(path)

        if self.sim:
            self.sim.reset(modname)
        else:
            self.sim = Simulator(modname, self.dt)
        #self.draw()
        self.paintGL()


    def animate(self):
        if self.sim:
            self.sim.step()
            self.updateSelectedCell()
            self.frameNo += 1
            if self.record:
                if (self.frameNo%5)==0:
                    #self.setSnapshotCounter(self.frameNo)
                    self.saveSnapshot()
    
    def saveSnapshot(self):
        print "saving snapshot", self.frameNo
        buf = self.grabFrameBuffer()
        return buf.save(os.path.join(
            self.sim.pickleDir,"{0:07}-frame.png".format(self.frameNo)), "png")

    def updateSelectedCell(self):
        if self.sim:
            states = self.sim.cellStates
            cid = self.selectedName
            txt = ''
            if states.has_key(cid):
                s = states[cid]
                for (name,val) in s.__dict__.items():
                    if name not in CellState.excludeAttr:
                        vals = str(val)
                        #if len(vals)>6: vals = vals[0:6]
                        txt = txt + name + ': ' + vals + '\n'
            self.selectedCell.emit(txt)
            if self.sim.stepNum%100==0:
                self.updateGL()

    def postSelection(self, name):
        self.selectedName = name
        self.updateSelectedCell()

    def paintGL(self):
        PyGLWidget.paintGL(self)
        glClearColor(0.5,0.5,0.5,0.0)
        glClear(GL_COLOR_BUFFER_BIT)
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        #s = self.renderInfo.scale
        #glScalef(s,s,s)
        if self.sim:
            for r in self.sim.renderers:
                if r != None:
                    r.render_gl(self.selectedName)

        glPopMatrix()

    def drawWithNames(self):
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        #s = self.renderInfo.scale
        #glScalef(s,s,s)
        if self.sim:
            for r in self.sim.renderers:
                if r:
                    r.renderNames_gl()
        glPopMatrix()


        
class RenderInfo:
    def __init__(self):
        self.renderers = []
        self.scale = 1.0
    def addRenderer(self, renderer):
        self.renderers.append(renderer)
    def reset(self):
        self.renderers = []
        self.scale = 1.0
    def setScale(self, s):
        self.scale = s

