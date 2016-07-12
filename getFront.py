#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
parser.add_argument("-v","--plaquevisibilitythreshold",type=float,default=1e6)
parser.add_argument("-f","--fractionthreshold",type=float,default=1e-4)
parser.add_argument("-z","--zerothreshold",type=float,default=1e-10)
args = parser.parse_args()

try:
    data = np.genfromtxt(args.infile)
except:
    print >> sys.stderr,"could not open file"
    exit(1)

timetrace = data[:,0]
timesteps = np.unique(timetrace)
 
# column:   2        3        4          5          6       7       8            9
methods = ['avginf','maxinf','avgphage','maxphage','bact1','bact2','visibility','susceptiblethreshold']

for t in timesteps:
    conf = data[timetrace == t]

    x = conf[:,1]
    nutrients = conf[:,2]
    b0 = conf[:,3]
    i0 = conf[:,4]
    b1 = conf[:,5]
    i1 = conf[:,6]
    p  = conf[:,7]
    
    b0[b0<args.zerothreshold] = 0
    b1[b1<args.zerothreshold] = 0
    i0[i0<args.zerothreshold] = 0
    i1[i1<args.zerothreshold] = 0
    p [p <args.zerothreshold] = 0
    
    allbact = b0 + i0 + b1 + i1
    s = np.zeros(len(b0))
    s[allbact>1] = b0[allbact>1] / allbact[allbact>1]

    pos = {}
    try:
        pos['avginf'] = np.dot(i0 + i1,x)/np.sum(i0 + i1)
    except:
        pos['avginf'] = 0
    try:
        pos['avgphage'] = np.dot(p,x)/np.sum(p)
    except:
        pos['avgphage'] = 0

    pos['maxinf']   = x[(i0+i1).argmax()]
    pos['maxphage'] = x[p.argmax()]
    
    pos['bact1']    = x[np.diff(b0,n=1).argmax()]
    pos['bact2']    = x[np.diff(b0,n=2).argmin()]
    
    allbactmt = allbact - args.plaquevisibilitythreshold
    if 0 < len(allbact[allbactmt < 0]) < len(allbactmt):
        idx1 = (allbactmt**2).argmin()
        if allbactmt[idx1] < 0:
            idxbelow = idx1
            if allbactmt[idx1+1] < 0:
                idxupper = idx1 - 1
            else:
                idxupper = idx1 + 1
        else:
            idxupper = idx1
            if allbactmt[idx1+1] < 0:
                idxbelow = idx1 + 1
            else:
                idxbelow = idx1 - 1
    try:
        pos['visibility'] = (x[idxupper] * allbactmt[idxupper] - x[idxbelow] * allbactmt[idxbelow])/(allbactmt[idxupper] - allbactmt[idxbelow])
    except:
        pos['visibility'] = 0
    
    pos['susceptiblethreshold'] = x[((s-args.fractionthreshold)**2).argmin()]
    #except:
        #pos['susceptiblethreshold'] = 0

    
    
    print "{:.3f}".format(t),
    for k in methods:
        print "{:.5f}".format(pos[k]),
    print

