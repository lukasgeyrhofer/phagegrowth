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
parser.add_argument("-x","--xpos",type=int,default=2850) # second lane = 2850, first lane = 2300, other parameters are identical
parser.add_argument("-X","--xwidth",type=int,default=100)
parser.add_argument("-y","--ypos",type=int,default=1000)
parser.add_argument("-Y","--ywidth",type=int,default=400)
parser.add_argument("-s","--smear",type=int,default=0,help="average profile over SMEAR pixels")
parser.add_argument("-m","--minpixel",type=int,default=380)
parser.add_argument("-M","--maxpixel",type=int,default=500)
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

    if args.smear > 1:
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
    #d,p = curve_fit(f,lum[i,:],lum[i+1,:],p0=np.ones(1))
    prof1 = lum[i,args.minpixel:args.maxpixel]
    prof2 = lum[i+1,args.minpixel:args.maxpixel]
    
    for diffcexp in xrange(-2,4,.1):
        print 10**diffcexp
    exit(1)
    
    interp1 = f(lum[i,:],1e-1)
    interp2 = f(lum[i,:],1e0)
    interp3 = f(lum[i,:],1e1)
    interp4 = f(lum[i,:],1e2)
    interp5 = f(lum[i,:],1e3)
    interp6 = f(lum[i,:],1e4)
    for j in range(len(lum[i,:])):
	print j,prof1[j],prof2[j],interp1[j],interp2[j],interp3[j],interp4[j],interp5[j],interp6[j]
    #print i,d[0]
    print





