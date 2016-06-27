#!/usr/bin/env python

import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
parser.add_argument("-c","--plaquevisibilitythreshold",type=float,default=1e6)
args = parser.parse_args()

try:
    data = np.genfromtxt(args.infile)
except:
    print >> sys.stderr,"could not open file"
    exit(1)

timetrace = data[:,0]
timesteps = np.unique(timetrace)

methods = ['avginf','maxinf','avgphage','maxphage','bact1','bact2','visibility']

for t in timesteps:
    conf = data[timetrace == t]

    x = conf[:,1]
    nutrients = conf[:,2]
    b0 = conf[:,3]
    i0 = conf[:,4]
    b1 = conf[:,5]
    i1 = conf[:,6]
    p  = conf[:,7]

    pos = {}
    pos['avginf']   = np.dot(i0 + i1,x)/np.sum(i0 + i1)
    pos['avgphage'] = np.dot(p,x)/np.sum(p)
    pos['maxinf']   = x[(i0+i1).argmax()]
    pos['maxphage'] = x[p.argmax()]
    
    pos['bact1']    = x[np.diff(b0,n=1).argmax()]
    pos['bact2']    = x[np.diff(b0,n=2).argmin()]
    
    allbact = b0 + i0 + b1 + i0 - args.plaquevisibilitythreshold
    if 0 < len(allbact[allbact < 0]) < len(allbact):
        idx1 = (allbact**2).argmin()
        if allbact[idx1] < 0:
            idxbelow = idx1
            if allbact[idx1+1] < 0:
                idxupper = idx1 - 1
            else:
                idxupper = idx1 + 1
        else:
            idxupper = idx1
            if allbact[idx1+1] < 0:
                idxbelow = idx1 + 1
            else:
                idxbelow = idx1 - 1
    try:
        pos['visibility'] = (x[idxupper] * allbact[idxupper] - x[idxbelow] * allbact[idxbelow])/(allbact[idxupper] - allbact[idxbelow])
    except:
        pos['visibility'] = 0
    
    
    print "{:.3f}".format(t),
    for k in methods:
        print "{:.5f}".format(pos[k]),
    print






