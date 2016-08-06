#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
from scipy.optimize import curve_fit


def printdata(x,y,stderr=False):
    if stderr:
        out = sys.stderr
    else:
        out = sys.stdout
    for d in zip(*(x,y)):
        print >> out,"{:6.4f} {:8.4f}".format(d[0],d[1])

parser = argparse.ArgumentParser()
parser.add_argument("-B","--BacterialGrowthInfile")
parser.add_argument("-P","--PhageGrowthInfile")
parser.add_argument("-m","--minDilution",default=0,type=float)
parser.add_argument("-M","--maxDilution",default=1,type=float)
parser.add_argument("-d","--dDilution",default=.01,type=float)
parser.add_argument("-X","--coefficients",default=[5.44661,26.1396,173.23],nargs="*",type=float)
args = parser.parse_args()


try:
    databg = np.genfromtxt(args.BacterialGrowthInfile)
    datapg = np.genfromtxt(args.PhageGrowthInfile)
except:
    print >> sys.stderr,"Could not open file"
    exit(1)

# estimate B-P growth curve from data
bactN   = databg[:,0]
bactGR  = databg[:,1]

phageN  = datapg[:,0]
phageGR = datapg[:,1]

Ephagebact = np.interp(bactN,phageN,phageGR)

# estimate B-P growth curve from quadratic fit for phage growth, P ~ Polynom(N^2,N,1)
# get polynomial fit
dilutionFIT            = np.arange(args.minDilution,args.maxDilution+args.dDilution,args.dDilution)
phageFIT               = np.zeros(len(dilutionFIT))
curN                   = np.ones(len(dilutionFIT))
PolynomialCoefficients = np.array(args.coefficients,dtype=float)

for coeff in PolynomialCoefficients:
    phageFIT += coeff * curN
    curN     *= dilutionFIT

EphagebactFIT = np.interp(bactN,dilutionFIT,phageFIT)

# in addition, interpolate growth curve of bacteria
EbactGR = np.interp(dilutionFIT,bactN,bactGR)




#printdata(bactGR,Ephagebact)
printdata(bactGR,EphagebactFIT,True)
printdata(EbactGR,phageFIT)

