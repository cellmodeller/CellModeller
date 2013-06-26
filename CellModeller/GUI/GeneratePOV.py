import sys
import os.path
import cPickle
import re
import math
from operator import itemgetter
import numpy

filename_regex = re.compile('-(\d+).pickle')

# For 3d
#preamble = """#version 3.6;
##declare RNG = seed(1234);
#global_settings{assumed_gamma 1.0}
#background{color rgb <0.0, 0.0, 0.0>}
#camera{location <0.0, -150.0, 100.0> direction z right x*image_width/image_height look_at <0.0, -0.0, 0.0>}
#light_source{<20,20,100> color rgb <0.8,0.8,0.8>}
#light_source{<0,0,200> color rgb <0.4,0.4,0.4> parallel point_at <0,0,0>}
#"""

# for xy slices
preamble = """#version 3.6;
declare RNG = seed(1234);
global_settings{assumed_gamma 1.0}
background{color rgb <0.0, 0.0, 0.0>}
camera{orthographic location <0,0,100> right <600,0,0> up <0,600,0> look_at <0,0,0>}
"""

# for xz slices
#preamble = """#version 3.6;
#declare RNG = seed(1234);
#global_settings{assumed_gamma 1.0}
#background{color rgb <0.0, 0.0, 0.0>}
#camera{orthographic location <0,100,0> right <200,0,0> up <0,0,25> look_at <0,0,0>}
#"""


capsule_str = """merge{\
cylinder{<%f,%f,%f>,<%f,%f,%f>,%f}\
sphere{<%f,%f,%f>,%f}\
sphere{<%f,%f,%f>,%f}\
texture{pigment{color rgb<%f,%f,%f>}}}\
"""

sliced_capsule_str = """intersection{\
merge{\
cylinder{<%f,%f,%f>,<%f,%f,%f>,%f}\
sphere{<%f,%f,%f>,%f}\
sphere{<%f,%f,%f>,%f}}\
plane{<%f,%f,%f>,0.01}\
plane{<-%f,-%f,-%f>,0.01}\
texture{pigment{color rgb<%f,%f,%f>}}\
}"""


def capsulestr(p, d, l, r, c):
    a = tuple(p - 0.5*l*d)
    b = tuple(p + 0.5*l*d)
    return capsule_str % (a + b + (r,) + a + (r,) + b + (r,) + c)


def get_color(species):
    r = 0
    g = 0.8
    b = 0.8
    return (r, g, b)

mxx = 0,0,0,0
mxy = 0,0,0,0
mxz = 0,0,0,0

def endscapsulestr(a,b,r,c):
    return capsule_str % (a + b + (r,) + a + (r,) + b + (r,) + c)

def sliced_capsulestr(pdir, pdist, a, b, r, c):
    return sliced_capsule_str % (a + b + (r,) + a + (r,) + b + (r,) + pdir + pdir + c)


def hollow_capsulestr(a, b, r, c, w):
    return hollow_capsule_str % (a + b + (r+w,) + a + (r+w,) + b + (r+w,) + a + b + (r,) + a + (r,) + b + (r,) + c)


def sliced_framestr(plane_dir, plane_dist, frameno, states):
    out = preamble
    discarded = 0
    for i,state in states.iteritems():
        a,b = state.ends
        a = tuple(a)
        b = tuple(b)
        r = state.radius
        #c = 0.0, 0.9, 0.9
        c = tuple(state.color)

        # for xz slice
        # ysameside = a[1]*b[1]
        # if ysameside > 0 and abs(a[1]) > r+0.2 and abs(b[1]) > r+0.2:
        #    discarded +=1
        #    continue

        # for xy slice
        if min(a[2], b[2]) > r+0.2:
            discarded += 1
            continue

        out += sliced_capsulestr(plane_dir, plane_dist, a, b, r-0.1, c)

    print 'discarded:', discarded
    return out


def framestr(frameno,states):
    out = preamble
    for i,state in states.iteritems():
        a,b = state.ends
        a = tuple(a)
        b = tuple(b)

        global mxx,mxy,mxz
        aa = a + (state.id,frameno, i, b, state.pos, state.dir, state.length, state.radius)
        bb = b + (state.id,frameno,)
        mxx = max(mxx, aa, bb, key=lambda x: abs(x[0]))
        mxy = max(mxy, aa, bb, key=lambda x: abs(x[1]))
        mxz = max(mxz, aa, bb, key=lambda x: abs(x[2]))

        r = state.radius
        #c = get_color(state.species)
        #c = 0.25,0.25,0.25
        c = tuple(state.color)
        out += endscapsulestr(a,b,r,c) + '\n'
    return out

def load_data(infile):
    return cPickle.load(infile)

def parse_filename(filename):
    print "filename = %s"%filename
    return int(filename_regex.search(filename).group(1))

def lineage(parents, founders, id):
    while id not in founders:
        id = parents[id]
    return id

def main(argv):
    indirname = argv[1]
    infilenames = os.listdir(indirname)
    infile_nums = [parse_filename(inf) for inf in infilenames]
    infilenames = map(lambda fn: os.path.join(indirname, fn), infilenames)
    infilenames = sorted(zip(infile_nums, infilenames))

    outdirname = argv[2]

    for frame,fn in infilenames:
        print 'reading file:',fn
        infile = open(fn, 'rb')
        states,parents = load_data(infile)
        infile.close()

        # color lineages
        for id,state in states.iteritems():
            p = lineage(parents, [1,2], id)
            if p == 1: state.color = [0, 1, 0]
            elif p == 2: state.color = [0, 0, 1]

        # for biophysCL-10-51-05-06-12/povxy
        outstr = sliced_framestr((0,0,1), 0, frame, states)

        # for biophysCL-10-51-05-06-12/povxz
        #outstr = sliced_framestr((0,1,0), 0, frame, states)

        #outstr = framestr(frame, states)
        outfile = open(os.path.join(outdirname,'frame%04i.pov'%(frame/10)), 'w')
        outfile.write(outstr)
        outfile.close()

    global mxx,mxy,mxz
    print mxx,mxy,mxz


if __name__ == '__main__':
    main(sys.argv)
