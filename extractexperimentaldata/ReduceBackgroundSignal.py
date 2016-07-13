#!/usr/bin/env python

import numpy as np
import argparse
import sys,math


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infiles",nargs="*")
parser.add_argument("-b","--backgroundpixel",default=100,type=int)
parser.add_argument("-o","--outsuffix",default=".reducedBG")
parser.add_argument("-m","--rescalewithmaximum",type=int,default=None)
args = parser.parse_args()


for fn in args.infiles:
    try:
        data = np.genfromtxt(fn)
    except:
        continue
    
    print fn
    p = data[:,0]
    l = data[:,1]
    
    plin = p[-args.backgroundpixel:]
    llin = l[-args.backgroundpixel:]
    space = len(plin)
    
    a = (space * np.sum(plin * llin) - np.sum(plin) * np.sum(llin)) / (space * np.sum(plin * plin) - np.sum(plin) * np.sum(plin))
    b = (np.sum(llin) - a * np.sum(plin))/space
    
    l_redbg = l - (a*p + b)
    
    if not args.rescalewithmaximum is None:
        l_redbg /= l_redbg[args.rescalewithmaximum]
    
    np.savetxt(fn + args.outsuffix,np.transpose([p,l_redbg]))
    
    







