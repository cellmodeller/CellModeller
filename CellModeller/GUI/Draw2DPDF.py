import sys
import math
import cPickle
import string
from reportlab.pdfgen import canvas
from reportlab.lib import units
import numpy

bg_color = (0.0,0.0,0.0)


def setup_canvas(name, world, page, page_center):
    c = canvas.Canvas(name)
    worldx,worldy = world
    pagex,pagey = page
    page_centx,page_centy = page_center
    c.setPageSize((pagex*units.cm, pagey*units.cm))
    c.setFillColorRGB(*bg_color)
    c.rect(0, 0, pagex*units.cm, pagey*units.cm, fill=1)
    c.translate(pagex*units.cm/2.0, pagey*units.cm/2.0)
    c.translate(page_centx*units.cm, page_centy*units.cm)
    c.scale(float(pagex*units.cm)/worldx, float(pagey*units.cm)/worldy)
    return c

def capsule_path(c, l, r):
    path = c.beginPath()
    path.moveTo(-l/2.0, -r)
    path.lineTo(l/2.0, -r)

    path.arcTo(l/2.0-r, -r, l/2.0+r, r, -90, 180)

    path.lineTo(-l/2.0, r)
    path.arc(-l/2.0-r, -r, -l/2.0+r, r, 90, 180)
    #path.close()
    return path

mx0 = 0
mx1 = 0
mx2 = 0
mx3 = 0
mx4 = 0
mnx = 10000
mxx = -10000
mny = 10000
mxy = -10000
def calc_cell_colors(state):
    global mx0,mx1,mx2,mx3,mx4,mnx,mxx,mny,mxy
    # global max_cfp, max_yfp
    # if state.cellType == 1:
    #     cfp = state.species[4]/8.
    #     max_cfp = max(max_cfp, cfp)
    #     r,g,b = 0,cfp,cfp
    # else:
    #     yfp = state.species[3]
    #     max_yfp = max(max_yfp, yfp)
    #     r,g,b = yfp,yfp,0
    # return (r,g,b),(0,0,0)
    # mx0 = max(mx0, state.species[0])
    # mx1 = max(mx1, state.species[1])
    # mx2 = max(mx2, state.species[2])
    # mx3 = max(mx3, state.species[3])
    # mx4 = max(mx4, state.species[4])
    # mnx = min(mnx, state.pos[0])
    # mxx = max(mxx, state.pos[0])
    # mny = min(mny, state.pos[1])
    # mxy = max(mxy, state.pos[1])



    # r = 0.2#0.1 + 0.9*state.species[0]/2300.0
    # g = 0.2 + 0.8*state.species[2]/0.4#0.1 + 0.9*state.species[1]/8000.0
    # b = 0.2 + 0.8*state.species[2]/0.4


    # # for EdgeDetectorChamber-22-57-02-06-12
    # # note i've reversed the colors from those in the file (makes no diff)
    # yfp = state.species[3]/3.34
    # cfp = state.species[4]/9.58
    # r = 0.2 + 0.8*yfp
    # g = 0.2 + 0.8*cfp + 0.8*yfp
    # b = 0.2 + 0.8*cfp


    # for Meristem-13-48-05-06-12
    #r = 0.2
    #g = 1.0
    #b = 0.2
    #if state.growthRate < 0.01:
    #    b = 1.0
    #else:
    #    r = 1.0


    '''colors = [ \
            [0.3, 0, 0], \
            [0, 1, 0], \
            [0, 0, 0.6], \
            [0.6, 0, 0], \
            [0.1, 0.3, 1], \
            [0, 0, 1], \
            [0, 0.6, 0], \
            [0, 0, 0.3], \
            [0.8, 0.1, 0.1] \
                ]'''
    colors = {0: (0.10255651959235379, 0.38540694651225138, 0.23256142392596002), 1: (0.28780219231559834, 0.3978802482183228, 0.032431247687302989), 2: (0.33929078662947099, 0.89102062747038346, 0.15162088685527164), 3: (0.85282408386397002, 0.067616659494312614, 0.54584642728543664), 4: (0.1613769193883654, 0.81478496354508712, 0.92119810819455483), 5: (0.45938783583329024, 0.14071040961916725, 0.71585065136216997), 6: (0.098193718762362603, 0.98468818000341485, 0.63713926021762912), 7: (0.60836498224956359, 0.16829106259571236, 0.15938712405805877), 8: (0.37101158138091073, 0.37644045628916378, 0.013125673969703433)}

    #colors = [ (0,1,0), (0,0,1) ]

    #(r,g,b) = state.color

    #if state.cellType!=0:
    #    r,g,b = colors[state.cellType]
    #    #r,g,b = state.color
    #else:
    #    r,g,b = (0,0.2,0.5)

    #seccolmap = {(0.2,0.2,0.2):(0.8,0.8,0.1), (0.2,0.8,0.2):(0.8,0.2,0.2), (0.2,0.2,0.8):(0.2,0.8,0.2)}
    #if tuple(state.color)==(0.2,0.2,0.2):
    #    return (0.2,0.8,0.8), (0.2,0.8,0.8)
    #else:
    #    r,g,b = state.color
    #(r,g,b) = seccolmap[tuple(state.color)]
    #return (r,g,b), (r,g,b)

    (r,g,b) = colors[state.cellType]
    return (r,g,b), (r,g,b)




def draw_capsule(c, p, d, l, r, fill_color, stroke_color):
    c.saveState()
    c.setStrokeColorRGB(*stroke_color)
    c.setFillColorRGB(*fill_color)
    c.setLineWidth(0.002*units.cm)
    c.setLineWidth(0.000*units.cm)
    c.translate(p[0], p[1])
    c.rotate(math.degrees(math.atan2(d[1], d[0])))
    path = capsule_path(c, l, r)
    c.drawPath(path, fill=1)
    c.restoreState()

def draw_cells(c, states):
    for id, state in states.items():
        p = state.pos
        d = state.dir
        l = state.length
        r = state.radius
        fill, stroke = calc_cell_colors(state)
        draw_capsule(c, p, d, l, r, fill, stroke)

def draw_chamber(c):
    # for EdgeDetectorChamber-22-57-02-06-12
    c.setLineWidth(0.015*units.cm)
    c.line(-100, -16, 100, -16)
    c.line(-100, 16, 100, 16)


mxsig0 = 0
def draw_signals(c, signals):
    global mxsig0
    # for EdgeDetectorChamber-22-57-02-06-12
    l, orig, dim, levels = signals
    l = map(float,l)
    for i in range(dim[1]):
        x = l[0]*i + orig[0]
        for j in range(dim[2]):
            y = l[1]*j + orig[1]
            lvl = levels[0][i][j][2]/0.0129
            mxsig0=max(lvl, mxsig0)
            c.setFillColorRGB(1.0-lvl, 1.0-lvl, 1.0-lvl)
            c.rect(x-l[0]/2.0, y-l[1]/2.0, l[0], l[1], stroke=0, fill=1)



def draw_frame(name, world, page, center, states, signals):
    c = setup_canvas(name, world, page, center)
    #draw_chamber(c)
    if signals: draw_signals(c, signals)
    draw_cells(c, states)
    c.showPage()
    c.save()

def lineage(parents, founders, id):
    while id not in founders:
        id = parents[id]
    return id


infns = sys.argv[1:]
for infn in infns:
    if infn[-4:]=='.pdf': continue
    frameno = int(infn[-11:-7])/10
    outfn = string.replace(infn, '.pickle', '-1um.pdf')
    print outfn
    data = cPickle.load(open(infn, 'rb'))
    if len(data) == 3:
        states, signals, parents = data
    else:
        states, parents = data
        signals = None

    #for id,state in states.iteritems():
    #    p = lineage(parents, [1,2], id)
    #    if p == 1: state.color = [0, 1, 0]
    #    elif p == 2: state.color = [1, 0, 0]

    # for EdgeDetectorChamber-22-57-02-06-12
    #world = (257.74,257.74)
    #world = (257.74/2.0,257.74/2.0)
    world = (550,550)
    #world = (710,710)
    #world = (155,155)
    page = (20,20)
    center = (0,0)

    # for Meristem-13-48-05-06-12
    #world = (424,62)
    #page = (21.2,3.1)
    #center = (-10.1,0)


    draw_frame(outfn, world, page, center, states, signals)



print mx0,mx1,mx2,mx3,mx4
print mnx, mxx, mny, mxy
print mxsig0
