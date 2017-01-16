#!/usr/bin/env python

import numpy as np
import argparse
import sys
from scipy.optimize import curve_fit

def radius(s,delta):
    herdimmunity = params['burstsize']*s - 1. - params['bacterialgrowthrate']*params['latentperiod']
    herdimmunity[herdimmunity < 0] = 0
    return 2*np.sqrt(params['diffusionconstant'] * delta * herdimmunity)*params['timetodepletion']


def radius2(s,delta,betaoverlambda):
    herdimmunity = betaoverlambda*s - 1
    herdimmunity[herdimmunity < 0] =0
    return 2*np.sqrt(params['diffusionconstant'] * delta * herdimmunity)*params['timetodepletion']
    

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
parser.add_argument("-D","--diffusionconstant",type=float,default=1.17e-2)
parser.add_argument("-b","--burstsize",type=float,default=85.6)
parser.add_argument("-l","--latentperiod",type=float,default=0.6017)
parser.add_argument("-a","--bacterialgrowthrate",type=float,default=.63)
parser.add_argument("-T","--timetodepletion",type=float,default=9.6)
parser.add_argument("-C","--excludecontrol",action="store_true",default=False)
parser.add_argument("-m","--maxfev",default=2000,type=int)
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
try:
    e = data[:,2]
except:
    pass

if args.excludecontrol:
    s = s[s<1]
    r = r[s<1]


initialguess = np.array([1e-2])
paramfit,paramcov = curve_fit(radius,s,r,p0 = initialguess,maxfev = args.maxfev)

print "# ************************************************************ #"
print "#   estimate absorption constant delta for phage on bacteria   #"
print "# ************************************************************ #"
print "# delta          = {:.6e}".format(paramfit[0])
print "# delta stddev   = {:.6e}".format(np.sqrt(paramcov[0][0]))
predictedradius = radius(s,paramfit[0])
print "# ************************************************************ #"
for key,value in params.iteritems():
    print "# {} = {}".format(key,value)
print "# radius(s) = 2*sqrt({:.6e} * {:.6e} * ({:.6e}*s - {:.6e} ))*{:.6e}".format(params['diffusionconstant'],paramfit[0],params['burstsize'],1-params['latentperiod']*params['bacterialgrowthrate'], params['timetodepletion'])
print "# ************************************************************ #"
print "#   data                                                       #"
print "# ************************************************************ #"
for i in range(len(s)):
    print "{:.2f} {:.6f} {:.6f}".format(s[i],predictedradius[i],r[i]),
    try:
        print "{:.6f}".format(e[i]),
    except:
        pass
    print


