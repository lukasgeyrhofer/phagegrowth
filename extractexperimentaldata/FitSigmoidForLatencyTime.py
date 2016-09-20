#!/usr/bin/env python

import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erf
import argparse
import sys,math


def sigmoid(x,center,slope,minval2,maxval2):
    return minval2**2 + (maxval2**2-minval2**2)/(1. + np.exp(-(x-center)/slope))

def errorfunction(x,center,slope,minval2,maxval2):
    return minval2**2 + (maxval2**2 - minval2**2) *0.5*(1+ erf((x-center)/(2.*slope)))

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
parser.add_argument("-m","--maxfev",type=int,default=5000)
args = parser.parse_args()

try:
    data = np.genfromtxt(args.infile)
except:
    print >> sys.stderr,"could not open file"
    exit(1)

t = data[:,0]
p = data[:,1]

initialguess = np.array([np.mean(t),10,p[0],p[-1]])

sfit,scov = curve_fit(sigmoid,t,p,p0=initialguess,maxfev=args.maxfev)
efit,ecov = curve_fit(errorfunction,t,p,p0=initialguess,maxfev=args.maxfev)
print "sigmoid:   {:.6f} {:.6f} {:.6e} {:.6e}".format(sfit[0],sfit[1],sfit[2]**2,sfit[3]**2)
print "errorfunc: {:.6f} {:.6f} {:.6e} {:.6e}".format(efit[0],efit[1],efit[2]**2,efit[3]**2)
print "sigm(t) = {} + {}/(1 + exp(-(t-{})/{}))".format(sfit[2]**2,sfit[3]**2-sfit[2]**2,sfit[0],sfit[1])
print "erfunc(t) = {} + {}*(1+erf((t-{})/{}))".format(efit[2]**2,.5*(efit[3]**2 - efit[2]**2),efit[0],2*efit[1])
