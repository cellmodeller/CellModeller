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
import cPickle

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
        self.set_radius(32)
        self.frameNo = 0
        self.fromP = False

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

    @pyqtSlot()
    def reset(self):
        if self.fromP==False:
            reload(self.sim.module)
        del self.sim

        if self.fromP==True:
            self.sim = Simulator(self.modname, self.dt, fromPickle=self.moduleStr)
        else:
            self.sim = Simulator(self.modname, self.dt)
        self.frameNo = 0

    @pyqtSlot()
    def loadPickle(self):
        qs = QtGui.QFileDialog.getOpenFileName(self, 'Load pickle file', '', '*.pickle')
        if qs:
            self.loadedPickle = str(qs)
            data = cPickle.load(open(self.loadedPickle,'rb'))
            if isinstance(data, dict):
                if self.sim:
                    del self.sim
                self.modname = data['moduleName']
                self.moduleStr = data['moduleStr']
                self.sim = Simulator(self.modname, self.dt, fromPickle=self.moduleStr)
                self.fromP = True
                self.sim.loadFromPickle(data)
                self.paintGL()
            else:
                print "Pickle is in an unsupported format, sorry"

    @pyqtSlot()
    def load(self):
        qs = QtGui.QFileDialog.getOpenFileName(self, 'Load Python module', '', '*.py')
        if qs:
            self.modfile = str(qs)
            self.loadFile(self.modfile)
            self.fromP = False

    def loadFile(self, modstr):
        (path,name) = os.path.split(modstr)
        modname = str(name).split('.')[0]
        self.modname = modname
        sys.path.append(path)
        if self.sim:
            del self.sim
        self.sim = Simulator(modname, self.dt)
        #self.draw()
        self.paintGL()


    def animate(self):
        if self.sim:
            self.sim.step()
            self.updateSelectedCell()
            self.frameNo += 1
    
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
            #if self.sim.stepNum%100==0:
            self.updateGL()

    def postSelection(self, name):
        self.selectedName = name
        self.updateSelectedCell()

    def translate(self, _trans):
        # Translate the selected cell by _trans (see PyGLWidget)
        cid = self.selectedName
        # Uncomment the condition below to enable movement of cells
        # currently the translation is not great due to the way PyGLWidget works
        if False: #self.sim and self.sim.cellStates.has_key(cid):
            self.sim.moveCell(cid, _trans)
            print "Called self.sim.moveCell"
            self.updateSelectedCell()
        else:
            # Translate the object by _trans
            # Update modelview_matrix_
            self.makeCurrent()
            glMatrixMode(GL_MODELVIEW)
            glLoadIdentity()
            glTranslated(_trans[0], _trans[1], _trans[2])
            glMultMatrixd(self.modelview_matrix_)
            self.modelview_matrix_ = glGetDoublev(GL_MODELVIEW_MATRIX)
            self.translate_vector_[0] = self.modelview_matrix_[3][0]
            self.translate_vector_[1] = self.modelview_matrix_[3][1]
            self.translate_vector_[2] = self.modelview_matrix_[3][2]
            self.signalGLMatrixChanged.emit()

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

