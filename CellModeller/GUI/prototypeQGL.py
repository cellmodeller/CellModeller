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
from PyQGLViewer import *

from OpenGL.GL import *
from OpenGL.GLU import *
import random
import os
import numpy
import time

import CellModeller.GUI.Renderers
from CellModeller import Simulator
from CellModeller.GUI.CMViewer import CMViewer, RenderInfo

import sys


qapp = QApplication([])
ui = uic.loadUi('CellModeller/GUI/test.ui')
ui.show()
cmv = ui.CMViewer
#cmv.setRegulator(reg)
#cmv.setRenderInfo(renderinfo)
sys.exit(qapp.exec_())
