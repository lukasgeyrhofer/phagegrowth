#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
args = parser.parse_args()

try:
    data = np.genfromtxt(args.infile)
except:
    print >> sys.stderr,"could not open file"
    exit(1)


idx   = data[:,0]
diffc = data[:,1]
dprof = data[:,2]

idxu  = np.unique(idx)

diffc = np.reshape(diffc,(len(idxu),len(idx)/len(idxu)))[0]
dprof = np.reshape(dprof,(len(idxu),len(idx)/len(idxu)))

#print diffc,idxu

for j in range(len(idxu)):
    print idxu[j],diffc[dprof[j].argmin()]


