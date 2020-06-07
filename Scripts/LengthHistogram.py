import sys
import os
import math
import numpy as np
sys.path.append('.')
import pickle
import CellModeller
#import matplotlib.pyplot as plt

file = sys.argv[1] #select pickle file in command line
output_file = open('LengthData.csv','w')


def rad_pos(cellstate):
    return np.sqrt(cellstate.pos[0]*cellstate.pos[0]+cellstate.pos[1]*cellstate.pos[1])

def lengthHist(pickle, bins, file=False):
    print(('opening '+ pickle))
    data = pickle.load(open(pickle,'r'))
    cs = data['cellStates']
    it = iter(cs)
    n = len(cs)
    print(('Number of cells = '+str(n)))
    lens = []
    r = []
    for it in cs:
        radialPosition = rad_pos(cs[it])
        cellLength = cs[it].length+2*cs[it].radius
        r.append(radialPosition)
        lens.append(cellLength)
        if file:
            file.write(str(radialPosition)+', '+ str(cellLength)+'\n')
    #plt.hist(lens,bins)
    #heatmap, xedges, yedges = np.histogram2d(r, lens, bins=100)
    #extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    #plt.imshow(heatmap,extent=extent,aspect='auto')
    #plt.show()
    #print r, lens
    return r, lens

def dirLengthHist(dir,bins,file=False):
    for f in os.listdir(dir):
        if '.pickle' in f:
            number = f[5:-7]
            if file:
                fout = open('LengthData'+number+'.csv')
            print(('step number = ', number))
            n, lens = lengthHist(dir + f,bins,file)
            if file:
                fout.close()

bin = np.arange(0.0,6.0,0.1)
lengthHist(file,bin,output_file)
output_file.close()

