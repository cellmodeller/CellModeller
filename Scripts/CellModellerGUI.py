#
# CellModeller4
#
# GUI
#
# Tim Rudge
# Jan 2011

from PyQt4.QtGui import QApplication
from PyQt4 import uic

import CellModeller.GUI.Renderers
from CellModeller import Simulator
from CellModeller.GUI.PyGLCMViewer import PyGLCMViewer, RenderInfo
from pkg_resources import resource_stream

import os
import sys

# The Qt application
qapp = QApplication([])

# The UI
uifile = resource_stream('CellModeller.GUI', 'PyGLGUI.ui')
ui = uic.loadUi(uifile)
ui.show()
ui.raise_()
cmv = ui.PyGLCMViewer

# Load a model if specified
if len(sys.argv) > 1: cmv.loadFile(sys.argv[1])

# Launch app main loop
sys.exit(qapp.exec_())
