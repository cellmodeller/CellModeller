#
# CellModeller4
#
# GUI
#
# Tim Rudge
# Jan 2011

from PyQt4.QtGui import QApplication
from PyQt4.QtCore import *
from PyQt4 import uic

from OpenGL.GL import *
from OpenGL.GLU import *
import random
import os
import numpy
import time

import CellModeller.GUI.Renderers
from CellModeller import Simulator
from CellModeller.GUI.PyGLCMViewer import PyGLCMViewer, RenderInfo

import sys


qapp = QApplication([])
ui = uic.loadUi('CellModeller/GUI/PyGLViewer.ui')
ui.show()
cmv = ui.PyGLCMViewer
#cmv.setRegulator(reg)
#cmv.setRenderInfo(renderinfo)
if len(sys.argv) > 1: cmv.loadFile(sys.argv[1])
sys.exit(qapp.exec_())
