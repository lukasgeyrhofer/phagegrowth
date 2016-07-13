#!/usr/bin/env python

import numpy as np
import argparse
import sys,math



def heatkernel(space,d):
    k = np.ones((space,space))/np.sqrt(4.*math.pi*d)
    for i in range(space):
	for j in range(space):
	    k[i,j] *= np.exp(-0.25*(i-j)**2/d)
    return k

def f(x,d):
    space = len(x)
    return np.dot(heatkernel(space,d),x)

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infiles",nargs="*")
parser.add_argument("-m","--minpixel",type=int,default=335)
parser.add_argument("-M","--maxpixel",type=int,default=595)
parser.add_argument("-b","--border",type=int,default=None) # final results: using -b 60
parser.add_argument("-f","--framejump",type=int,default=1)
parser.add_argument("-d","--minlog10diffc",type=float,default=-3)
parser.add_argument("-D","--maxlog10diffc",type=float,default=2)
parser.add_argument("-S","--steplog10diffc",type=float,default=.05)
args = parser.parse_args()


if not len(args.infiles) > 1:
    print >> sys.stderr,"not enough files"
    exit(1)
else:
    data = np.genfromtxt(args.infiles[0])
    lastprofile = data[args.minpixel:args.maxpixel,1]
    
    i = args.framejump
    while i < len(args.infiles):
        data = np.genfromtxt(args.infiles[i])
        curprofile = data[args.minpixel:args.maxpixel,1]
        
        #compprofile01 = f(lastprofile,11.2014*5)
        #compprofile02 = f(lastprofile,0.10)
        #compprofile03 = f(lastprofile,0.12)
        #compprofile04 = f(lastprofile,0.14)
        #compprofile05 = f(lastprofile,0.16)
        #compprofile06 = f(lastprofile,0.20)
        #compprofile07 = f(lastprofile,0.50)
        #compprofile08 = f(lastprofile,1.00)
        #compprofile09 = f(lastprofile,2.00)
        #compprofile10 = f(lastprofile,5.00)
        #compprofile11 = f(lastprofile,10.00)
        #compprofile12 = f(lastprofile,20.00)
        #compprofile13 = f(lastprofile,50.00)
        #compprofile14 = f(lastprofile,100.00)
        #compprofile15 = f(lastprofile,200.00)
        #compprofile16 = f(lastprofile,500.00)
        #compprofile17 = f(lastprofile,1000.00)
        
        #for j in range(len(curprofile)):
            #print j,curprofile[j],lastprofile[j],compprofile01[j] #,compprofile02[j],compprofile03[j],compprofile04[j],compprofile05[j],compprofile06[j],compprofile07[j],compprofile08[j],compprofile09[j],compprofile10[j],compprofile11[j],compprofile12[j],compprofile13[j],compprofile14[j],compprofile15[j],compprofile16[j],compprofile17[j]
        #exit(1)
        
        basedifference = np.sum((lastprofile[args.border:-args.border] - curprofile[args.border:-args.border])**2)
        for diffconstexp in np.arange(args.minlog10diffc,args.maxlog10diffc,args.steplog10diffc):
            print "{} {:.6e} {:.10e}".format(i,10**diffconstexp,np.sum( (f(lastprofile,10**diffconstexp)[args.border:-args.border] - curprofile[args.border:-args.border])**2)/basedifference)
        
        lastprofile = curprofile[:]
        i+=args.framejump
        print
        
        
        
