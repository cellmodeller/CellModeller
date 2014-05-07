# -*- coding: utf-8 -*-
#===============================================================================
#
# PyGLWidget.py
#
# A simple GL Viewer.
#
# Copyright (c) 2011, Arne Schmitz <arne.schmitz@gmx.net>
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the <organization> nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#===============================================================================

from PyQt4 import QtCore, QtGui, QtOpenGL
import math
import numpy
import numpy.linalg as linalg
import OpenGL
OpenGL.ERROR_CHECKING = True
from OpenGL.GL import *
from OpenGL.GLU import *

class PyGLWidget(QtOpenGL.QGLWidget):

    # Qt signals
    signalGLMatrixChanged = QtCore.pyqtSignal()
    rotationBeginEvent = QtCore.pyqtSignal()
    rotationEndEvent = QtCore.pyqtSignal()

    def __init__(self, parent = None):
        format = QtOpenGL.QGLFormat()
        format.setSampleBuffers(True)
        QtOpenGL.QGLWidget.__init__(self, format, parent)
        #self.setCursor(QtCore.Qt.OpenHandCursor)
        self.setMouseTracking(True)

        self.modelview_matrix_  = []
        self.translate_vector_  = [0.0, 0.0, 0.0]
        self.viewport_matrix_   = []
        self.projection_matrix_ = []
        self.near_   = 0.1
        self.far_    = 100.0
        self.fovy_   = 45.0
        self.radius_ = 5.0
        self.last_point_2D_ = QtCore.QPoint()
        self.last_point_ok_ = False
        self.last_point_3D_ = [1.0, 0.0, 0.0]
        self.isInRotation_  = False
        self.pickSize = 3

        # connections
        #self.signalGLMatrixChanged.connect(self.printModelViewMatrix)

    @QtCore.pyqtSlot()
    def printModelViewMatrix(self):
        print self.modelview_matrix_

    def initializeGL(self):
        # OpenGL state
        glClearColor(0.0, 0.0, 0.0, 0.0)
        glEnable(GL_DEPTH_TEST)
        self.reset_view()

    def resizeGL(self, width, height):
        glViewport( 0, 0, width, height );
        self.set_projection( self.near_, self.far_, self.fovy_ );
        self.updateGL()

    def paintGL(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        glMatrixMode(GL_MODELVIEW)
        glLoadMatrixd(self.modelview_matrix_)

    def set_projection(self, _near, _far, _fovy):
        self.near_ = _near
        self.far_ = _far
        self.fovy_ = _fovy
        self.makeCurrent()
        glMatrixMode( GL_PROJECTION )
        glLoadIdentity()
        gluPerspective( self.fovy_, float(self.width()) / float(self.height()),
                        self.near_, self.far_ )
        self.updateGL()

    def set_pick_projection(self, x, y, _near, _far, _fovy):
        self.near_ = _near
        self.far_ = _far
        self.fovy_ = _fovy
        self.makeCurrent()
        glMatrixMode( GL_PROJECTION )
        glLoadIdentity()
        viewport =glGetIntegerv(GL_VIEWPORT)
        gluPickMatrix(x, viewport[3]-y, self.pickSize, self.pickSize, viewport);
        gluPerspective( self.fovy_, float(self.width()) / float(self.height()),
                        self.near_, self.far_ )
    
    def set_center(self, _cog):
        self.center_ = _cog
        self.view_all()

    def set_radius(self, _radius):
        self.radius_ = _radius
        self.set_projection(_radius / 100.0, _radius * 100.0, self.fovy_)
        self.reset_view()
        self.translate([0, 0, -_radius * 2.0])
        self.view_all()
        self.updateGL()

    def reset_view(self):
        # scene pos and size
        glMatrixMode( GL_MODELVIEW )
        glLoadIdentity();
        self.modelview_matrix_ = glGetDoublev( GL_MODELVIEW_MATRIX )
        self.set_center([0.0, 0.0, 0.0])

    def reset_rotation(self):
        self.modelview_matrix_[0] = [1.0, 0.0, 0.0, 0.0]
        self.modelview_matrix_[1] = [0.0, 1.0, 0.0, 0.0]
        self.modelview_matrix_[2] = [0.0, 0.0, 1.0, 0.0]
        glMatrixMode(GL_MODELVIEW)
        glLoadMatrixd(self.modelview_matrix_)
        self.updateGL()
   
    def translate(self, _trans):
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

    def rotate(self, _axis, _angle):
        t = [self.modelview_matrix_[0][0] * self.center_[0] +
             self.modelview_matrix_[1][0] * self.center_[1] +
             self.modelview_matrix_[2][0] * self.center_[2] +
             self.modelview_matrix_[3][0],
             self.modelview_matrix_[0][1] * self.center_[0] +
             self.modelview_matrix_[1][1] * self.center_[1] +
             self.modelview_matrix_[2][1] * self.center_[2] +
             self.modelview_matrix_[3][1],
             self.modelview_matrix_[0][2] * self.center_[0] +
             self.modelview_matrix_[1][2] * self.center_[1] +
             self.modelview_matrix_[2][2] * self.center_[2] +
             self.modelview_matrix_[3][2]]

        self.makeCurrent()
        glLoadIdentity()
        glTranslatef(t[0], t[1], t[2])
        glRotated(_angle, _axis[0], _axis[1], _axis[2])
        glTranslatef(-t[0], -t[1], -t[2])
        glMultMatrixd(self.modelview_matrix_)
        self.modelview_matrix_ = glGetDoublev(GL_MODELVIEW_MATRIX)
        self.signalGLMatrixChanged.emit()

    def view_all(self):
        self.translate( [ -( self.modelview_matrix_[0][0] * self.center_[0] +
                             self.modelview_matrix_[0][1] * self.center_[1] +
                             self.modelview_matrix_[0][2] * self.center_[2] +
                             self.modelview_matrix_[0][3]),
                           -( self.modelview_matrix_[1][0] * self.center_[0] +
                              self.modelview_matrix_[1][1] * self.center_[1] +
                              self.modelview_matrix_[1][2] * self.center_[2] +
                              self.modelview_matrix_[1][3]),
                           -( self.modelview_matrix_[2][0] * self.center_[0] +
                              self.modelview_matrix_[2][1] * self.center_[1] +
                              self.modelview_matrix_[2][2] * self.center_[2] +
                              self.modelview_matrix_[2][3] +
                              self.radius_ / 2.0 )])

    def map_to_sphere(self, _v2D):
        _v3D = [0.0, 0.0, 0.0]
        # inside Widget?
        if (( _v2D.x() >= 0 ) and ( _v2D.x() <= self.width() ) and
            ( _v2D.y() >= 0 ) and ( _v2D.y() <= self.height() ) ):
            # map Qt Coordinates to the centered unit square [-0.5..0.5]x[-0.5..0.5]
            x  = float( _v2D.x() - 0.5 * self.width())  / self.width()
            y  = float( 0.5 * self.height() - _v2D.y()) / self.height()

            _v3D[0] = x;
            _v3D[1] = y;
            # use Pythagoras to comp z-coord (the sphere has radius sqrt(2.0*0.5*0.5))
            z2 = 2.0*0.5*0.5-x*x-y*y;
            # numerical robust sqrt
            _v3D[2] = math.sqrt(max( z2, 0.0 ))

            # normalize direction to unit sphere
            n = linalg.norm(_v3D)
            _v3D = numpy.array(_v3D) / n

            return True, _v3D
        else:
            return False, _v3D

    def wheelEvent(self, _event):
        # Use the mouse wheel to zoom in/out
        
        d = - float(_event.delta()) / 200.0 * self.radius_
        self.translate([0.0, 0.0, d])
        self.updateGL()
        _event.accept()

    def selectName(self, point):
        glSelectBuffer(50) # allocate a selection buffer of SIZE elements
        glRenderMode(GL_SELECT)
        
        glMatrixMode( GL_PROJECTION )
        glPushMatrix()
        self.set_pick_projection( point.x(), point.y(), self.near_, self.far_, self.fovy_ );
        
        self.paintGL()
        self.drawWithNames()

        buf = glRenderMode(GL_RENDER)
        selectedName = -1
        for hit_record in buf:
            min_depth, max_depth, names = hit_record
            for name in names:
                if name:
                    selectedName = name
        glMatrixMode( GL_PROJECTION )
        glPopMatrix()
        return selectedName
    
    def mousePressEvent(self, _event):
        self.last_point_2D_ = _event.pos()
        self.last_point_ok_, self.last_point_3D_ = self.map_to_sphere(self.last_point_2D_)
        if (_event.buttons() & QtCore.Qt.LeftButton) and (_event.modifiers() & QtCore.Qt.ShiftModifier):
            name = self.selectName(_event.pos())
            self.postSelection(name)

    def mouseMoveEvent(self, _event):
        newPoint2D = _event.pos()

        if ((newPoint2D.x() < 0) or (newPoint2D.x() > self.width()) or
            (newPoint2D.y() < 0) or (newPoint2D.y() > self.height())):
            return
        
        # Left button: rotate around center_
        # Middle button: translate object
        # Left & middle button: zoom in/out

        value_y = 0
        newPoint_hitSphere, newPoint3D = self.map_to_sphere(newPoint2D)

        dx = float(newPoint2D.x() - self.last_point_2D_.x())
        dy = float(newPoint2D.y() - self.last_point_2D_.y())

        w  = float(self.width())
        h  = float(self.height())

        # enable GL context
        self.makeCurrent()

        # move in z direction
        if (((_event.buttons() & QtCore.Qt.LeftButton) and (_event.buttons() & QtCore.Qt.MidButton))
            or (_event.buttons() & QtCore.Qt.LeftButton and _event.modifiers() & QtCore.Qt.ControlModifier)):
            value_y = self.radius_ * dy * 2.0 / h;
            self.translate([0.0, 0.0, value_y])
        # move in x,y direction
        elif (_event.buttons() & QtCore.Qt.RightButton):
            z = - (self.modelview_matrix_[0][2] * self.center_[0] +
                   self.modelview_matrix_[1][2] * self.center_[1] +
                   self.modelview_matrix_[2][2] * self.center_[2] +
                   self.modelview_matrix_[3][2]) / (self.modelview_matrix_[0][3] * self.center_[0] +
                                                    self.modelview_matrix_[1][3] * self.center_[1] +
                                                    self.modelview_matrix_[2][3] * self.center_[2] +
                                                    self.modelview_matrix_[3][3])

            fovy   = 45.0
            aspect = w / h
            n      = 0.01 * self.radius_
            up     = math.tan(fovy / 2.0 * math.pi / 180.0) * n
            right  = aspect * up

            self.translate( [2.0 * dx / w * right / n * z,
                             -2.0 * dy / h * up / n * z,
                             0.0] )

    
        # rotate
        elif (_event.buttons() & QtCore.Qt.LeftButton and (not _event.modifiers() & QtCore.Qt.ShiftModifier)):
            if (not self.isInRotation_):
                self.isInRotation_ = True
                self.rotationBeginEvent.emit()
       
            axis = [0.0, 0.0, 0.0]
            angle = 0.0

            if (self.last_point_ok_ and newPoint_hitSphere):
                axis = numpy.cross(self.last_point_3D_, newPoint3D)
                cos_angle = numpy.dot(self.last_point_3D_, newPoint3D)
                if (abs(cos_angle) < 1.0):
                    angle = math.acos(cos_angle) * 180.0 / math.pi
                    angle *= 2.0
                self.rotate(axis, angle)

        # remember this point
        self.last_point_2D_ = newPoint2D
        self.last_point_3D_ = newPoint3D
        self.last_point_ok_ = newPoint_hitSphere

        # trigger redraw
        self.updateGL()

        def mouseReleaseEvent(self, _event):
            if (isInRotation_):
                isInRotation_ = false
                self.rotationEndEvent.emit()
            last_point_ok_ = False

#===============================================================================
#
# Local Variables:
# mode: Python
# indent-tabs-mode: nil
# End:
#
#===============================================================================
