import sys
import os
import math
import numpy as np
sys.path.append('.')
import pickle
import CellModeller
import matplotlib.pyplot as plt

pname= sys.argv[1]
print(('opening '+ pname))
fout=open('spatial.txt', "w")
fout.write("r spec_level \n")

bin_num=20
rad_max=95

narray=np.zeros(bin_num)
specArray=np.zeros(bin_num)

data = pickle.load(open(pname,'r'))
cs = data['cellStates']
it = iter(cs)
n = len(cs)

rArray = np.multiply(list(range(0,bin_num)),rad_max/bin_num)

for it in cs:
    r = np.sqrt(cs[it].pos[0]*cs[it].pos[0]+cs[it].pos[1]*cs[it].pos[1])
    for x in range(0,bin_num):
        if((r>=x*rad_max/bin_num)and(r<(x+1)*rad_max/bin_num)):
            specArray[x]=specArray[x]*narray[x]+cs[it].species[0]
            narray[x]+=1
            specArray[x]=specArray[x]/narray[x]

for x in range(0,bin_num):
    fout.write(str((x+1)*rad_max/bin_num/2)+" "+str(narray[x]))
    if(narray[x]==0): fout.write(" 0.0\n")
    else: fout.write(" "+str(specArray[x])+"\n")

fout.close()

plt.plot(rArray,specArray)
plt.show()