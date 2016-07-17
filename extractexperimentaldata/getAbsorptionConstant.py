#!/usr/bin/env python

import numpy as np
import argparse
import sys
from scipy.optimize import curve_fit

def radius(s,delta):
    return 2*np.sqrt(params['diffusionconstant'] * delta * (params['burstsize']/params['latentperiod']*s - 1))/params['bacterialgrowthrate']*np.log(params['bacterialgrowthratio'])


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
parser.add_argument("-D","--diffusionconstant",type=float,default=1.17e-2)
parser.add_argument("-b","--burstsize",type=float,default=91)
parser.add_argument("-l","--latentperiod",type=float,default=.5)
parser.add_argument("-a","--bacterialgrowthrate",type=float,default=.72)
parser.add_argument("-B","--bacterialgrowthratio",type=float,default=1e3)
parser.add_argument("-C","--excludecontrol",action="store_true",default=False)
args = parser.parse_args()



try:
    data = np.genfromtxt(args.infile)
except:
    print >> sys.stderr,"could not open file"
    exit(1)

global params
params = vars(args)

s = data[:,0]
r = data[:,1]

if args.excludecontrol:
    s = s[s<1]
    r = r[s<1]


deltafit,deltacov = curve_fit(radius,s,r)

print "# ************************************************************ #"
print "#   estimate absorption constant delta for phage on bacteria   #"
print "# ************************************************************ #"
print "# delta = {:.10e}".format(deltafit[0])
print "# deltastddev = {:.10e}".format(np.sqrt(deltacov[0][0]))
print "# ************************************************************ #"
for key,value in params.iteritems():
    print "# {} = {}".format(key,value)
print "# diameter(s) = 4*sqrt({:.10e} * {:.10e} * (({:.10e}/{:.10e} + 1)*s -1 ))/{:.10e}*log({:.10e})".format(params['diffusionconstant'],deltafit[0],params['burstsize'],params['latentperiod'],params['bacterialgrowthrate'],params['bacterialgrowthratio'])
print "# ************************************************************ #"
print "#   data                                                       #"
print "# ************************************************************ #"
for i in range(len(s)):
    print "{:.2f} {:.10f} {:.10f}".format(s[i],radius(s[i],deltafit[0]),r[i])


