#!/usr/bin/env python

import numpy as np
import argparse
import sys
from scipy.optimize import curve_fit

def radius(s,delta):
    herdimmunity = params['burstsizeoverlatentperiod']*s - 1
    herdimmunity[herdimmunity < 0] =0
    return 2*np.sqrt(params['diffusionconstant'] * delta * herdimmunity)/params['bacterialgrowthrate']*np.log(params['bacterialgrowthratio'])


def radius2(s,delta,betaoverlambda):
    herdimmunity = betaoverlambda*s - 1
    herdimmunity[herdimmunity < 0] =0
    return 2*np.sqrt(params['diffusionconstant'] * delta * herdimmunity)/params['bacterialgrowthrate']*np.log(params['bacterialgrowthratio'])
    

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
parser.add_argument("-D","--diffusionconstant",type=float,default=1.17e-2)
parserBL = parser.add_mutually_exclusive_group()
parserBL.add_argument("-E","--estimatebetaoverlambda",default=False,action="store_true")
parserBL.add_argument("-b","--burstsizeoverlatentperiod",type=float,default=182.)
parser.add_argument("-a","--bacterialgrowthrate",type=float,default=.72)
parser.add_argument("-B","--bacterialgrowthratio",type=float,default=1e3)
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


if args.estimatebetaoverlambda:
    initialguess = np.array([1e-2,2e2])
    paramfit,paramcov = curve_fit(radius2,s,r,p0 = initialguess,maxfev = args.maxfev)
else:
    initialguess = np.array([1e-2])
    paramfit,paramcov = curve_fit(radius,s,r,p0 = initialguess,maxfev = args.maxfev)

print "# ************************************************************ #"
print "#   estimate absorption constant delta for phage on bacteria   #"
print "# ************************************************************ #"
print "# delta          = {:.10e}".format(paramfit[0])
print "# delta stddev   = {:.10e}".format(np.sqrt(paramcov[0][0]))
if args.estimatebetaoverlambda:
    print "# betaoverlambda = {:.10e}".format(paramfit[1])
    print "# bl stddev      = {:.10e}".format(np.sqrt(paramcov[1][1]))
    params['burstsizeoverlatentperiod'] = paramfit[1]
    predictedradius = radius2(s,paramfit[0],paramfit[1])
else:
    predictedradius = radius(s,paramfit[0])
print "# ************************************************************ #"
for key,value in params.iteritems():
    print "# {} = {}".format(key,value)
print "# radius(s) = 2*sqrt({:.10e} * {:.10e} * ({:.10e}*s -1 ))/{:.10e}*log({:.10e})".format(params['diffusionconstant'],paramfit[0],params['burstsizeoverlatentperiod'],params['bacterialgrowthrate'],params['bacterialgrowthratio'])
print "# ************************************************************ #"
print "#   data                                                       #"
print "# ************************************************************ #"
for i in range(len(s)):
    print "{:.2f} {:.10f} {:.10f}".format(s[i],predictedradius[i],r[i]),
    try:
        print "{:.10f}".format(e[i]),
    except:
        pass
    print


