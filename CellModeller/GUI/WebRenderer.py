"""
I'm putting this in its own file to avoid importing all the OpenGL stuff 
"""
import numpy
import math

class WebRenderer:
    class Sphere:
        def __init__(self, pos, rad, color):
            assert type(pos) is list
            assert len(pos) == 3
            
            assert type(color) is list
            assert len(color) == 3 or len(color) == 4
        
            self.position = pos
            self.radius = rad
            self.color = color
            
    class SignalsGrid:
        def __init__(self, origin, cell_size, cell_count, colors_data):
            self.origin = [ origin[i] for i in range(3) ]
            self.cell_size = [ cell_size[i] for i in range(3) ]
            self.cell_count = [ cell_count[i] for i in range(3) ]
            self.voxels = colors_data

    def __init__(self):
        self.spheres = []
            
        self.sig = None
        self.integ = None
        self.has_signals_grid = False
    
    def addSphere(self, pos, radius, color=[ 0.8, 0.8, 0.8, 0.8 ]):
        self.spheres.append(self.Sphere(pos, radius, color))
        
    def attachSignals(self, sig, integ):
        self.sig = sig
        self.integ = integ
        self.has_signals_grid = True
    
    def getSignalsGrid(self):
        if not self.has_signals_grid:
            return None
        
        dim = self.sig.gridDim[1:]
        
        signalLevel = self.integ.signalLevel.reshape(self.sig.gridDim)
        imageData = signalLevel[0]
        
        mx = numpy.max(imageData)
        mn = numpy.min(imageData)
        scale = (255 / (mx-mn) ) if mx > mn else 1
        
        imageData = (imageData - mn) * scale
        
        return self.SignalsGrid(self.sig.gridOrig, self.sig.gridSize, dim, imageData)