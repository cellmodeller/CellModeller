import PyQt4
from PyQt4 import QtCore, QtGui
from PyQt4.Qt import Qt
from PyQt4.QtCore import QObject, QTimer, pyqtSignal, pyqtSlot, QStringList
from PyGLWidget import PyGLWidget
from OpenGL.GL import *
from OpenGL.GLU import *

from CellModeller.Regulation import ModuleRegulator
from CellModeller.Simulator import Simulator
from CellModeller.CellState import CellState
import os
import sys
import cPickle
import pyopencl as cl

class PyGLCMViewer(PyGLWidget):

    setSavePicklesToggle = pyqtSignal(bool)
    selectedCell = pyqtSignal(str) #emit selected cell info
    selectedName = -1
    dt = 0.05

    def __init__(self, parent = None):
        PyGLWidget.__init__(self,parent)
        self.animTimer = QTimer()
        self.animTimer.timeout.connect(self.animate)
        self.animTimer.start(0)
        self.renderInfo = None
        self.sim = None
        self.run = False
        self.frameNo = 0
        self.loadingFromPickle = False
        self.clPlatformNum=0
        self.clDeviceNum=0

        # Initial view setup
        self.set_radius(32)
        self.rotate([0,0,1],-45)
        self.translate([0,0,20])
        self.rotate([1,0,0],-45)

    def help(self):
        pass

    def setSimulator(self, sim):
        if self.sim:
            del self.sim
        self.sim = sim
        self.frameNo = sim.stepNum
        if self.run:
            self.frameNo += 1
        # Make GUI button match simulator state for saving pickles
        self.setSavePicklesToggle.emit(sim.saveOutput)
        # Get rid of any selected cell id
        self.selectedName = -1

    def getOpenCLPlatDev(self):
        return self.getOpenCLPlatform() and self.getOpenCLDevice()

    def getOpenCLPlatform(self):
        # Pop dialogs to get user to choose OpenCL platform 
        platforms = cl.get_platforms()

        platlist = QStringList([str(p.name) for p in platforms])
        platdict = dict(zip(platlist, range(len(platlist))))

        if len(platlist)==1:
            self.clPlatformNum = 0
            return True

        qsPlatformName, ok = QtGui.QInputDialog.getItem(self, \
                                            'Choose OpenCL platform', \
                                            'Available platforms:', \
                                            platlist, \
                                            editable=False)
        if not ok:
            print "You didn't select a OpenCL platform..."
            return False
        else:
            self.clPlatformNum = platdict[qsPlatformName]
            return True

    def getOpenCLDevice(self):
        # Pop dialogs to get user to choose OpenCL platform and device
        platforms = cl.get_platforms()
        devices = platforms[self.clPlatformNum].get_devices()

        devlist = QStringList([str(d.name) for d in devices])
        devdict = dict(zip(devlist, range(len(devlist))))
        
        if len(devlist)==1:
            self.clDeviceNum = 0
            return True
        
        qsDeviceName, ok = QtGui.QInputDialog.getItem(self, \
                                            'Choose OpenCL device', \
                                            'Available devices:', \
                                            devlist, \
                                            editable=False)
        if not ok:
            print "You didn't select a OpenCL device..."
            return False
        else:
            self.clDeviceNum = devdict[qsDeviceName]
            return True


    @pyqtSlot(bool)
    def toggleRun(self, run):
        self.run = run
        if run:
            self.frameNo += 1

    @pyqtSlot(bool)
    def toggleSavePickles(self, save):
        self.writePickles = save
        self.sim.setSaveOutput(save)

    @pyqtSlot()
    def reset(self):
        # Note: we don't ask user to choose OpenCL platform/device on reset, only load
        if not self.loadingFromPickle:
            reload(self.sim.module)

        if self.loadingFromPickle:
            sim = Simulator(self.modName, \
                            self.dt, \
                            moduleStr=self.moduleStr, \
                            clPlatformNum=self.clPlatformNum, \
                            clDeviceNum=self.clDeviceNum, \
                            is_gui=True) 
            self.setSimulator(sim) 
        else:
            sim = Simulator(self.modName, \
                                self.dt, \
                                clPlatformNum=self.clPlatformNum, \
                                clDeviceNum=self.clDeviceNum, \
                                is_gui=True) 
            self.setSimulator(sim) 
        self.frameNo = 0
        self.updateGL()

    @pyqtSlot()
    def loadPickle(self):
        qs = QtGui.QFileDialog.getOpenFileName(self, 'Load pickle file', '', '*.pickle')
        if qs and self.getOpenCLPlatDev():
            filename = str(qs)
            data = cPickle.load(open(filename,'rb'))
            if isinstance(data, dict):
                self.modName = data['moduleName']
                self.moduleStr = data['moduleStr']
                self.frameNo = data['stepNum']
                sim = Simulator(self.modName, \
                                    self.dt, \
                                    moduleStr=self.moduleStr, \
                                    clPlatformNum=self.clPlatformNum, \
                                    clDeviceNum=self.clDeviceNum, \
                                    is_gui=True) 
 
                self.setSimulator(sim) 
                self.loadingFromPickle = True
                self.sim.loadFromPickle(data)
                # Note: the pickle loaded contains the stepNum, hence we now
                # need to set the GUI frameNo to match
                self.frameNo = self.sim.stepNum
                if self.run:
                    self.frameNo += 1
                self.updateGL()
            else:
                print "Pickle is in an unsupported format, sorry"

    @pyqtSlot()
    def load(self):
        qs = QtGui.QFileDialog.getOpenFileName(self, 'Load Python module', '', '*.py')
        if qs:
            modfile = str(qs)
            self.loadModelFile(modfile)

    def loadModelFile(self, modname):
        if self.getOpenCLPlatDev():
            self.loadingFromPickle=False
            (path,name) = os.path.split(modname)
            self.modName = str(name).split('.')[0]
            sys.path.append(path)
            sim = Simulator(self.modName, \
                                    self.dt, \
                                    clPlatformNum=self.clPlatformNum, \
                                    clDeviceNum=self.clDeviceNum, \
                                    is_gui=True) 
 
            self.setSimulator(sim)
            self.updateGL()


    def animate(self):
        if self.sim:
            if self.frameNo > self.sim.stepNum:
                if self.sim.step():
                    self.updateSelectedCell()
                    if self.run:
                        self.frameNo += 1
    
    def updateSelectedCell(self):
        if self.sim:
            states = self.sim.cellStates
            cid = self.selectedName
            txt = ''
            if states.has_key(cid):
                txt += 'Selected Cell (id = %d)\n---\n'%(cid)
                s = states[cid]
                for (name,val) in s.__dict__.items():
                    if name not in CellState.excludeAttr:
                        txt += name + ': '
                        txt += str(val)
                        txt += '\n'
            self.selectedCell.emit(txt)
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

        # Draw a grid in xy plane
        glDisable(GL_LIGHTING)
        glColor3f(1.0, 1.0, 1.0)
        glEnable(GL_LINE_SMOOTH)
        glLineWidth(1.0)
        glBegin(GL_LINES)
        for i in range(5):
            glVertex(-20, (i-2)*10)
            glVertex(20, (i-2)*10)
            glVertex((i-2)*10, -20)
            glVertex((i-2)*10, 20)
        glEnd()

        # Draw x,y,z axes
        glDisable(GL_DEPTH_TEST)
        glBegin(GL_LINES)
        glColor3f(1.0,0.0,0.0)
        glVertex(0,0,0)
        glVertex(5,0,0)
        glColor3f(0.0,1.0,0.0)
        glVertex(0,0,0)
        glVertex(0,5,0)
        glColor3f(0.0,0.0,1.0)
        glVertex(0,0,0)
        glVertex(0,0,5)
        glEnd()
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_LIGHTING)

        # Draw model
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

