#!/usr/bin/env python

from PIL import Image
import numpy as np
import argparse
import sys,math

from scipy.optimize import curve_fit


def heatkernel(space,d):
    k = np.ones((space,space))/np.sqrt(4*math.pi*d)
    for i in range(space):
	for j in range(space):
	    k[i,j] *= np.exp(-0.25*(i-j)**2/d)
    return k

def f(x,d):
    space = len(x)
    return np.dot(heatkernel(space,d),x)


def gray(r,g,b):
    # luminosity, computes grayscale by different weighting to adjust for human perception
    return (0.21*r + 0.72*g + 0.07*b)/256.

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile",nargs="*")
parser.add_argument("-x","--xpos",type=int,default=2300)
parser.add_argument("-X","--xwidth",type=int,default=100)
parser.add_argument("-y","--ypos",type=int,default=1000)
parser.add_argument("-Y","--ywidth",type=int,default=400)
parser.add_argument("-s","--smear",type=int,default=0,help="average profile over SMEAR pixels")
args = parser.parse_args()

imgcount = len(args.infile)
if imgcount <= 1:
    print >> sys.stderr,"can not compute time evolution, not enough images!"
    exit(1)
lum      = np.zeros((imgcount,args.ywidth))
lumsmear = np.zeros((imgcount,args.ywidth))


# read pictures
for i in range(imgcount):
    print >> sys.stderr,"# reading file '%s'"%args.infile[i]
    img = Image.open(args.infile[i])

    lum2d = np.zeros((args.xwidth,args.ywidth))
    lum1d = np.zeros(args.ywidth)

    for y in range(args.ywidth):
        for x in range(args.xwidth):
            r,g,b = img.getpixel((x+args.xpos,y+args.ypos))
            lum2d[x,y] = gray(r,g,b)
        lum1d[y] = np.mean(lum2d[:,y])
        lum[i,:] = lum1d[:]

    if args.smear > 0:
	for y in range(args.ywidth):
	    smear1 = max(0,y-int(args.smear/2))
	    smear2 = min(y+int(args.smear/2),args.ywidth)
	    lumsmear[i,y] = np.mean(lum1d[smear1:smear2])
	lum[i,:] = lumsmear[i,:]
    

# code to evaluate diffusion equation directly

#expandlum = np.concatenate((np.reshape(lumsmear[:,0],newshape=(imgcount,1)),lumsmear,np.reshape(lumsmear[:,-1],newshape=(imgcount,1))),axis=1)
#d2xlumtmp = np.diff(expandlum,n=2,axis=1)
#d2xlum = np.zeros((imgcount-1,args.ywidth))
#for i in range(imgcount-1):
    #d2xlum[i] = 0.5*d2xlumtmp[i] + 0.5*d2xlumtmp[i+1]
#dtlum = np.diff(lumsmear,axis=0)


# code to propagate solutions with heat kernel and estimate diffusion constant from this propagation
for i in range(imgcount-1):
    print >> sys.stderr,"# analyzing transition (%d)->(%d)"%(i,i+1)
    d,p = curve_fit(f,lum[i,:],lum[i+1,:],p0=np.ones(1))
    print i,d[0]
    





