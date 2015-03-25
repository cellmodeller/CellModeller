import sys
import math
import cPickle
import string
from reportlab.pdfgen.canvas import Canvas
from reportlab.lib import units
import numpy
from reportlab.lib.colors import Color
import math


class CellModellerPDFGenerator(Canvas):
    # ---
    # Class that extends reportlab pdf canvas to draw CellModeller simulations
    # ---
    def __init__(self, name, data):
        self.name = name
        self.states = data.get('cellStates')
        self.signals = data.get('signals')
        self.parents = data.get('lineage')
        self.data = data
        Canvas.__init__(self, name)

    def setup_canvas(self, name, world, page, page_center, bg_color):
        #c = canvas.Canvas(name)
        worldx,worldy = world
        pagex,pagey = page
        page_centx,page_centy = page_center
        self.setPageSize((pagex*units.cm, pagey*units.cm))
        self.setFillColor(bg_color)
        self.rect(0, 0, pagex*units.cm, pagey*units.cm, fill=1)
        self.translate(pagex*units.cm/2.0, pagey*units.cm/2.0)
        self.translate(page_centx*units.cm, page_centy*units.cm)
        self.scale(float(pagex*units.cm)/worldx, float(pagey*units.cm)/worldy)
        #return c

    def capsule_path(self, l, r):
        path = self.beginPath()
        path.moveTo(-l/2.0, -r)
        path.lineTo(l/2.0, -r)

        path.arcTo(l/2.0-r, -r, l/2.0+r, r, -90, 180)

        path.lineTo(-l/2.0, r)
        path.arc(-l/2.0-r, -r, -l/2.0+r, r, 90, 180)
        #path.close()
        return path

    def calc_cell_colors(self, state):
        # Generate Color object from cellState, with black outline
        (r,g,b) = state.color
        return [Color(r,g,b,alpha=1.0) , Color(0,0,0,alpha=1.0)]

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
        for id, state in self.states.items():
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

    mxsig0 = 0
    def draw_signals(self):
        global mxsig0
        # for EdgeDetectorChamber-22-57-02-06-12
        l, orig, dim, levels = self.signals
        l = map(float,l)
        for i in range(dim[1]):
            x = l[0]*i + orig[0]
            for j in range(dim[2]):
                y = l[1]*j + orig[1]
                lvl = levels[0][i][j][2]/0.0129
                mxsig0=max(lvl, mxsig0)
                self.setFillColorRGB(1.0-lvl, 1.0-lvl, 1.0-lvl)
                self.rect(x-l[0]/2.0, y-l[1]/2.0, l[0], l[1], stroke=0, fill=1)

    def draw_frame(self, name, world, page, center, bg_color):
        self.setup_canvas(name, world, page, center, bg_color)
        #draw_chamber(c)
        if self.signals: 
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
        for (id,s) in self.states.iteritems():
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
        print 'Importing CellModeller pickle file: %s'%fname
        data = cPickle.load(open(fname, 'rb'))
        return data
    else:
        return None


def main():
    # To do: parse some useful arguments as rendering options here
    # e.g. outline color, page size, etself.
    #
    # For now, put these options into variables here:
    bg_color = Color(1.0,1.0,1.0,alpha=1.0)

    # For now just assume a list of files
    infns = sys.argv[1:]
    for infn in infns:
        # File names
        if infn[-4:]=='.pdf': continue
        frameno = int(infn[-11:-7])/10
        outfn = string.replace(infn, '.pickle', '-1um.pdf')
        print 'Frame %d: Processing %s to generate %s'%(frameno,infn,outfn)
        
        # Import data
        data = importPickle(infn)
        if not data:
            print "Problem importing data!"
            return

        # Create a pdf canvas thing
        pdf = CellModellerPDFGenerator(outfn, data)

        # Get the bounding square of the colony to size the image
        # This will resize the image to fit the page...
        # ** alternatively you can specify a fixed world size here
        (w,h) = pdf.computeBox()
        world = (w,h)

        # Page setup
        page = (20,20)
        center = (0,0)

        # Render pdf
        print 'Rendering PDF output to %s'%outfn
        pdf.draw_frame(outfn, world, page, center, bg_color)

if __name__ == "__main__": 
    main()

