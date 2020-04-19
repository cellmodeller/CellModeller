#
# CellModeller4
#
# GUI
#
# Tim Rudge
# Jan 2011

from PyQt5.QtWidgets import QApplication
from PyQt5 import uic
from PyQt5.QtCore import *

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
pix_ratio = qapp.devicePixelRatio()
cmv.setPixelRatio(pix_ratio)
label = ui.label
label.setTextFormat(Qt.RichText)
label.setAlignment(Qt.AlignJustify)

# Load a model if specified
if len(sys.argv) > 1: cmv.loadModelFile(sys.argv[1])

# Launch app main loop
sys.exit(qapp.exec_())
