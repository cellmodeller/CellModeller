import PyQt4
from PyQt4 import QtCore, QtGui
from PyQt4.Qt import Qt
from PyQt4.QtCore import QObject, pyqtSignal, pyqtSlot
from PyQGLViewer import *
from OpenGL.GL import *
from OpenGL.GLU import *

from CellModeller.Regulation import ModuleRegulator
from CellModeller.Simulator import Simulator
from CellModeller.CellState import CellState
import os
import sys


class CMViewer(QGLViewer):

    selectedCell = pyqtSignal(str)#CellState, name='selectedCell')

    def __init__(self, parent = None):
	QGLViewer.__init__(self,parent)
	hw = self.helpWidget()
	if hw:
	    hw.hide()
        self.renderInfo = None
	self.sim= None
	self.record = False
        self.setSceneRadius(1000)
        self.setSelectRegionWidth(1)
        self.setSelectRegionHeight(1)
        self.setAnimationPeriod(50)
        self.showEntireScene()
	self.setSnapshotFormat('JPEG')
	self.setSnapshotFileName('Images/test')
	self.setMouseBinding(Qt.LeftButton + Qt.RightButton, QGLViewer.CAMERA, QGLViewer.ZOOM)
	self.frameNo = 0
    
    def help(self):
	pass

    def setSimulator(self, sim):
	self.sim = sim

    @pyqtSlot(bool)
    def toggleRecord(self, rec):
	self.record = rec

    @pyqtSlot()
    def reset(self):
	if self.sim:
	    self.sim.reset()
	self.frameNo = 0
	
    @pyqtSlot()
    def load(self):
	qs =QtGui.QFileDialog.getOpenFileName(self, 'Load Python module', '', '*.py')
	modstr = str(qs)
	(path,name) = os.path.split(modstr)
	modname = str(name).split('.')[0]
	sys.path.append(path)

	if self.sim:
	    self.sim.reset(modname)
	else:
	    self.sim = Simulator(modname)
	self.draw()
	

    def animate(self):
	if self.sim:
	    self.sim.step(0.25)
	    self.updateSelectedCell()
	    self.frameNo += 1
	    if self.record:
		if (self.frameNo%5)==0:
		    self.setSnapshotCounter(self.frameNo)
		    self.saveSnapshot()

    def updateSelectedCell(self):
	states = self.sim.cellStates
	cid = self.selectedName()
	txt = ''
	if states.has_key(cid):
	    s = states[cid]
	    for (name,val) in s.__dict__.items():
		if name not in CellState.excludeAttr:
		    vals = str(val)
		    #if len(vals)>6: vals = vals[0:6]
		    txt = txt + name + ': ' + vals + '\n'

	self.selectedCell.emit(txt)
	self.updateGL()

    def postSelection(self, point):
	self.updateSelectedCell()

    def draw(self):
        self.clipPlane()

        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        #s = self.renderInfo.scale
        #glScalef(s,s,s)
	if self.sim:
            for r in self.sim.renderers:
                if r != None:
                    r.render_gl(self.selectedName())

        glPopMatrix()

    def drawWithNames(self):
        self.clipPlane()
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        #s = self.renderInfo.scale
        #glScalef(s,s,s)
        for r in self.sim.renderers:
            if r != None:
                r.renderNames_gl()

        glPopMatrix()

    def clipPlane(self):
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glMultMatrixd(self.manipulatedFrame().matrix())
        glClipPlane(GL_CLIP_PLANE0, [ 0.0, 0.0, 1.0, 0.0 ])
        glColor3f(0,0,0)
    #glScalef(1000.0,1000.0,1000.0)
        self.drawArrow(2, .5)
        glBegin(GL_LINE_STRIP)
        glVertex3f(-3.0, -3.0, 0.001)
        glVertex3f(-3.0,  3.0, 0.001)
        glVertex3f( 3.0,  3.0, 0.001)
        glVertex3f( 3.0, -3.0, 0.001)
        glVertex3f(-3.0, -3.0, 0.001)
        glEnd()
        glPopMatrix()
        
    def init(self):
        self.restoreStateFromFile()
        #self.help()
        self.setManipulatedFrame(ManipulatedFrame())
        glEnable(GL_CLIP_PLANE0)
	if self.sim:
	    for r in self.sim.renderers:
		if r != None:
		    r.init_gl()
        glClearColor(0.5,0.5,0.5,1)
    
    def closeEvent(self,event):
        helpwidget = self.helpWidget()
        if not helpwidget is None and helpwidget.isVisible() :
            helpwidget.hide()
        QGLViewer.closeEvent(self,event)
        
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

