#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import math,sys


def growthrate(nutrients):
    #   function that reduces cell's growth rate with depletion of nutrients
    #   Kmax/((1+Ks)/c_nutrients)
    return bacteria['growth_Kmax']*nutrients/(nutrients+bacteria['growth_Kc'])

def burstsize(nutrients):
    #   linear function that reduces phage burst size on susceptibles with respect to bacterial growth rate
    return (phage['burstsize_max']-phage['burstsize_min'])*growthrate(nutrients)/bacteria['growth_Kmax']+phage['burstsize_min']


def f(y,yd):
  # unpack arrays
  n,s,si,r,ri,p = y
  siDELAY,riDELAY = yd
  # assume phage concentration is contant at boundaries
  # shifted arrays necessary for diffusion term of phage movement
  pfwd  = np.concatenate((p[1:],np.array([p[-1]])))
  pback = np.concatenate((np.array([p[0]]),p[:-1]))
  
  gr = growthrate(n)
  bs = burstsize(n)

  return np.array([	-gr*bacteria['invyield']*(s+r), \
			gr*s - phage['absorption']*s*p, \
			phage['absorption'] *s*p - phage['burstrate']*siDELAY, \
			gr*r - phage['absorption']*bacteria['resistance']*p*r, \
			phage['absorption']*bacteria['resistance']*r*p - phage['burstrate'] * riDELAY ,\
			bs * phage['burstrate'] * (siDELAY+riDELAY) - phage['absorption']*(s+si+r+ri)*p + phage['diffusion']*(pback - 2*p + pfwd)])

def RungeKutta4(y,yd):
  # 4th order Runge-Kutta integration scheme
  #print f(y,yd[0])
  #exit(1)
  k1 = param['epsilon']*f(y,yd[0])
  k2 = param['epsilon']*f(y+k1/2.,(yd[0]+yd[1])/2.)
  k3 = param['epsilon']*f(y+k2/2.,(yd[0]+yd[1])/2.)
  k4 = param['epsilon']*f(y+k3,yd[1])
  ret = y + (k1+2*k2+2*k3+k4)/6.

  ret[ret<0] = 0
  return ret

def main():
    parser = argparse.ArgumentParser()
    
    # different argument groups allow additional information to be printed to the screen
    # in the end everything is still parsed to the same "args" variable
    
    # define lattice
    parser_lattice = parser.add_argument_group(description="Lattice parameters")
    parser_lattice.add_argument("-s","--space",type=int,default=250)        # number of bins
    parser_lattice.add_argument("-d","--dx",type=float,default=0.02)        # in [mm]
    
    # algorithm parameters
    parser_algorithm = parser.add_argument_group(description="Algorithm parameters")
    parser_algorithm.add_argument("-e","--epsilon",type=float,default=5e-3) # measured in hrs, thus default eps=5e-3h=0.3min
    parser_algorithm.add_argument("-T","--maxTime",type=int,default=10)     # in hrs
    parser_algorithm.add_argument("-O","--outputsteps",type=int,default=20) # 20 5e-3 = 0.1h = 6min
    parser_algorithm.add_argument("-o","--outputbins",type=int,default=1)   
    parser_algorithm.add_argument("-Q","--quiet",default=False,action="store_true")

    # interactions between various trophic levels
    parser_interactions = parser.add_argument_group(description="Interactions between different trophic levels")
    parser_interactions.add_argument("-D","--phage_diffusionconstant",type=float,default=2.5e-2)
    parser_interactions.add_argument("-b","--phage_burstdelay",type=float,default=.5)
    parser_interactions.add_argument("-w","--phage_delaydistr",type=float,default=.05)
    parser_interactions.add_argument("-L","--phage_burstsize_linear",type=float,default=100)
    parser_interactions.add_argument("-l","--phage_burstsize_const",type=float,default=2)
    parser_interactions.add_argument("-A","--phage_bacteria_absorptionconstant",type=float,default=1e-5)
    parser_interactions.add_argument("-m","--bacteria_growth_Kmax",type=float,default=0.7204)
    parser_interactions.add_argument("-c","--bacteria_growth_Kc",type=float,default=0.0000257024)
    parser_interactions.add_argument("-n","--nutrients_per_cell",type=float,default=1e-10)          # mg/cell
    parser_interactions.add_argument("-r","--bacteria_resistant_absorptionreduction",type=float,default=1e-4)

    # startconcentrations 
    parser_startconc = parser.add_argument_group(description="Start concentrations")
    parser_startconc.add_argument("-S","--startconcentrationsusceptibles",type=float,default=10e4) # Indiv/mm
    parser_startconc.add_argument("-R","--startconcentrationresistant",   type=float,default=0e4)  # Indiv/mm
    parser_startconc.add_argument("-P","--startconcentrationphages",      type=float,default=100)  # Indiv in first bin
    parser_startconc.add_argument("-N","--startconcentrationnutrients",   type=float,default=1e-4) # mg/mm
    parser_startconc.add_argument("-p","--initialplaquesize",             type=float,default=.5)   # radius in mm


    args = parser.parse_args()
    
    global bacteria
    bacteria = {'growth_Kmax' : args.bacteria_growth_Kmax, 'growth_Kc' : args.bacteria_growth_Kc, 'invyield' : args.nutrients_per_cell, 'resistance': args.bacteria_resistant_absorptionreduction}
    global phage
    phage = {'diffusion' : args.phage_diffusionconstant/args.dx**2, 'burstrate':1./args.phage_burstdelay, 'bursttime_mean' : args.phage_burstdelay, 'bursttime_stddev' : args.phage_delaydistr, 'burstsize_max' : args.phage_burstsize_linear, 'burstsize_min': args.phage_burstsize_const, 'absorption': args.phage_bacteria_absorptionconstant}
    global param
    param = {'epsilon' : args.epsilon, 'dx' : args.dx}
    
    print bacteria
    print phage
    print param
    #exit(1)
    
    # array with all concentrations/numbers of cells
    y = np.zeros((6,args.space))

    # store solutions to implement time-delay in phage burst
    delaysize = int((args.phage_burstdelay+3*args.phage_delaydistr)/args.epsilon)
    delaytime = np.arange(delaysize)*args.epsilon-3*args.phage_delaydistr
    # two distributions shifted by a timestep for use in Runge-Kutta algorithm
    delaydistr1 = np.exp(-0.5*(delaytime/args.phage_delaydistr)**2)/np.sqrt(2*math.pi*args.phage_delaydistr**2)*args.epsilon
    delaydistr2 = np.exp(-0.5*((delaytime+args.epsilon)/args.phage_delaydistr)**2)/np.sqrt(2*math.pi*args.phage_delaydistr**2)*args.epsilon
    
    sdelay = np.zeros((delaysize,args.space))
    rdelay = np.zeros((delaysize,args.space))
    
    y[0,:] = args.startconcentrationnutrients
    y[1,:] = args.startconcentrationsusceptibles
    y[2,:] = 0
    y[3,:] = args.startconcentrationresistant
    y[4,:] = 0
    # not concentrations, but rather numbers of cells/phages in a given bin
    y[5,:int(args.initialplaquesize/args.dx)] = args.startconcentrationphages
    
    maxsteps = int(args.maxTime/args.epsilon)

    for o in range(maxsteps):
	# delay as convolution with 
	yd = np.array([[np.dot(delaydistr1,sdelay),np.dot(delaydistr1,rdelay)],[np.dot(delaydistr2,sdelay),np.dot(delaydistr2,rdelay)]])
	
	# integrate system of differential equations with Runge-Kutta method
	y = RungeKutta4( y ,yd)
 	
	# update stored solutions:
	# drop first item (0 as python index), append new solutions from current timestep to the end
	sdelay[0:delaysize-1,:] = sdelay[1:delaysize,:]
	sdelay[delaysize-1,:] = y[2]
	rdelay[0:delaysize-1,:] = rdelay[1:delaysize,:]
	rdelay[delaysize-1,:] = y[4]

	if o%args.outputsteps == 0:
	    # integrate profiles to get total number of ...
	    nn = np.sum(y[0])*args.dx          # nutrients
	    ss = np.sum(y[1])*args.dx          # susceptibles
	    si = np.sum(y[2])*args.dx          # infected
	    rr = np.sum(y[3])*args.dx          # resistant
	    ri = np.sum(y[4])*args.dx          # res inf
	    pp = np.sum(y[5])*args.dx          # phages
	    print "%10.4lf %.10e %.10e %.10e %.10e %.10e %.10e"%(o*args.epsilon,nn,ss,si,rr,ri,pp)
	    
	    if not args.quiet:
		for i in range(0,args.space,args.outputbins):
		    print >> sys.stderr,"%lf %lf %.10e %.10e %.10e %.10e %.10e %.10e"%(o*args.epsilon,i*args.dx,y[0,i],y[1,i],y[2,i],y[3,i],y[4,i],y[5,i])
		print >> sys.stderr


# more "pythonian" version of the program:
# whole file can be loaded as module, and then be executed as function call to "main()" (as defined above)
# if executed on its own, then python automatically goes here (and runs "main()")
if __name__=="__main__":
  main()
  
  
  

  
