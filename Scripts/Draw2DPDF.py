import os
import sys
import math
import string

from reportlab.pdfgen.canvas import Canvas
from reportlab.lib import units
from reportlab.lib.colors import Color

import numpy
import pickle

mxsig0 = 0

class CellModellerPDFGenerator(Canvas):
    # ---
    # Class that extends reportlab pdf canvas to draw CellModeller simulations
    # ---
    def __init__(self, name, data, bg_color):
        self.name = name
        self.states = data.get('cellStates')

        self.signal_levels = data.get('sigGrid', None)
        self.signal_grid_orig = data.get('sigGridOrig', None)
        self.signal_grid_dim = data.get('sigGridDim', None)
        self.signal_grid_size = data.get('sigGridSize', None)
        self.signals = self.signal_levels is not None
        
        self.parents = data.get('lineage')
        self.data = data
        self.bg_color = bg_color
        Canvas.__init__(self, name)

    # ----
    # Inherit this class and override the following method to change
    # how cell color is computed from current CellState.
    # Default behaviour is to use CellState.color, and outline in black
    # ----
    def calc_cell_colors(self, state):
        # Generate Color objects from cellState, with black outline
        (r,g,b) = state.color
        # Return value is tuple of colors, (fill, stroke)
        return [Color(r,g,b,alpha=1.0) , Color(0,0,0,alpha=1.0)]
    # ----

    def setup_canvas(self, name, world, page, page_center):
        worldx,worldy = world
        pagex,pagey = page
        page_centx,page_centy = page_center
        self.setPageSize((pagex*units.cm, pagey*units.cm))
        self.setFillColor(self.bg_color)
        self.rect(0, 0, pagex*units.cm, pagey*units.cm, fill=1)
        self.translate(pagex*units.cm/2.0, pagey*units.cm/2.0)
        self.translate(page_centx*units.cm, page_centy*units.cm)
        self.scale(float(pagex*units.cm)/worldx, float(pagey*units.cm)/worldy)

    def capsule_path(self, l, r):
        path = self.beginPath()
        path.moveTo(-l/2.0, -r)
        path.lineTo(l/2.0, -r)

        path.arcTo(l/2.0-r, -r, l/2.0+r, r, -90, 180)

        path.lineTo(-l/2.0, r)
        path.arc(-l/2.0-r, -r, -l/2.0+r, r, 90, 180)
        #path.close()
        return path

    def draw_capsule(self, p, d, l, r, fill_color, stroke_color):
        self.saveState()
        self.setStrokeColor(stroke_color)
        self.setFillColor(fill_color)
        self.setLineWidth(0.001*units.cm)
        self.setLineWidth(0.001*units.cm)
        self.translate(p[0], p[1])
        self.rotate(math.degrees(math.atan2(d[1], d[0])))
        path = self.capsule_path(l, r)
        self.drawPath(path, fill=1)
        self.restoreState()

    def draw_cells(self):
        for id, state in list(self.states.items()):
            p = state.pos
            d = state.dir
            l = state.length
            r = state.radius
            fill, stroke = self.calc_cell_colors(state)
            self.draw_capsule(p, d, l, r, fill, stroke)

    def draw_chamber(self):
        # for EdgeDetectorChamber-22-57-02-06-12
        self.setLineWidth(0.015*units.cm)
        self.line(-100, -16, 100, -16)
        self.line(-100, 16, 100, 16)

    def draw_signals(self, index=0, scale=0.0192, z=0):
        '''
        Function for drawing signal grids, currently limited to 1 signal a plane at a fixed z-axis level through the
        grid

        index = index of signal to render
        scale = scale factor for signal level
        z = height of slice through grid 
        '''
        global mxsig0
        # for EdgeDetectorChamber-22-57-02-06-12
        l, orig, dim, levels = self.signal_grid_size, \
                                self.signal_grid_orig, \
                                self.signal_grid_dim, \
                                self.signal_levels
        levels = levels.reshape(dim)
        mx = levels[index,:,:,z].max()
        l = list(map(float,l))
        for i in range(dim[1]):
            x = l[0]*i + orig[0]
            for j in range(dim[2]):
                y = l[1]*j + orig[1]
                lvls = levels[index,i,j,z]/mx
                mxsig0 = max(lvls, mxsig0)
                self.setFillColorRGB(lvls, 0, 0)
                self.rect(x-l[0]/2.0, y-l[1]/2.0, l[0], l[1], stroke=0, fill=1)

    def draw_frame(self, name, world, page, center):
        self.setup_canvas(name, world, page, center)
        #draw_chamber(c)
        if self.signals: 
            print("Drawing signals")
            self.draw_signals()
        self.draw_cells()
        self.showPage()
        self.save()

    def lineage(self, parents, founders, id):
        while id not in founders:
            id = parents[id]
        return id

    def computeBox(self):
        # Find bounding box of colony, minimum size = (40,40)
        mxx = 20
        mnx = -20
        mxy = 20
        mny = -20
        for (id,s) in list(self.states.items()):
            pos = s.pos    
            l = s.length    # add/sub length to keep cell in frame
            mxx = max(mxx,pos[0]+l) 
            mnx = min(mnx,pos[0]-l) 
            mxy = max(mxy,pos[1]+l) 
            mny = min(mny,pos[1]-l) 
        w = (mxx-mnx)
        h = (mxy-mny)
        return (w,h)
#
#End class CellModellerPDFGenerator

def importPickle(fname):
    if fname[-7:]=='.pickle':
        print(('Importing CellModeller pickle file: %s'%fname))
        data = pickle.load(open(fname, 'rb'))

        # Check for old-style pickle that is tuple, 
        # just extract cellStates from 1st element
        if isinstance(data, tuple):
            data = {'cellStates':data[0]}
        # Return dictionary of simulation data
        return data
    else:
        return None

# Define a pdf generator class with cell outline color same as fill color
class MyPDFGenerator(CellModellerPDFGenerator):
    def calc_cell_colors(self, state):
        # Generate Color objects from cellState, fill=stroke
        (r,g,b) = state.color
        # Return value is tuple of colors, (fill, stroke)
        fcol = Color(r,g,b,alpha=1.0)
        scol = Color(r*0.5,g*0.5,b*0.5,alpha=1.0)
        return [fcol,scol]

def main():
    # To do: parse some useful arguments as rendering options here
    # e.g. outline color, page size, etc.
    #
    # For now, put these options into variables here:
    bg_color = Color(0,0,0,alpha=1.0)

    # For now just assume a list of files
    infns = sys.argv[1:]
    for infn in infns:
        # File names
        if infn[-7:]!='.pickle':
            print(('Ignoring file %s, because its not a pickle...'%(infn)))
            continue

        outfn = infn.replace('.pickle', '.pdf')
        outfn = os.path.basename(outfn) # Put output in this dir
        print(('Processing %s to generate %s'%(infn,outfn)))
        
        # Import data
        data = importPickle(infn)
        if not data:
            print("Problem importing data!")
            return

        # Create a pdf canvas thing
        pdf = MyPDFGenerator(outfn, data, bg_color)

        # Get the bounding square of the colony to size the image
        # This will resize the image to fit the page...
        # ** alternatively you can specify a fixed world size here
        '''(w,h) = pdf.computeBox()
        sqrt2 = math.sqrt(2)
        world = (w/sqrt2,h/sqrt2)'''
        world = (150,150)

        # Page setup
        page = (20,20)
        center = (0,0)

        # Render pdf
        print(('Rendering PDF output to %s'%outfn))
        pdf.draw_frame(outfn, world, page, center)

if __name__ == "__main__": 
    main()

