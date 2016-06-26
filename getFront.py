#!/usr/bin/env python

import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
args = parser.parse_args()

try:
    data = np.genfromtxt(args.infile)
except:
    print >> sys.stderr,"could not open file"
    exit(1)

timetrace = data[:,0]
timesteps = np.unique(timetrace)

methods = ['avginf','maxinf','avgphage','maxphage']

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
    
    print "{:.3f}".format(t),
    for k in methods:
        print "{:.5f}".format(pos[k]),
    print






