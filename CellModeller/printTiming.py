import sys
import os
sys.path.append('.')
import cPickle
import CellModeller

dir = os.path.join('data', sys.argv[1])

print "No. cells\ttime (s)"
print "-------------------"
for f in os.listdir(dir):
    if 'pickle' in f:
        ff = os.path.join(dir,f)
        (cs,lin) = cPickle.load(open(ff,'r'))
        n = len(cs)
        t = os.path.getmtime(ff)
        print "%i\t%f"%(n,t)

