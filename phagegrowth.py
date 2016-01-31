#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import math,sys


def growthrate(nutrients):
    #   function that reduces cell's growth rate with depletion of nutrients
    #   Kmax/((1+Ks)/c_nutrients)
    return args.bacteria_growth_Kmax*nutrients/(nutrients+args.bacteria_growth_Kc)
    # np.array([args.bacteria_growth_Kmax/(1+(args.bacteria_growth_Kc/n)) if n>0 else 0 for n in nutrients])

def burstsize(nutrients):
    #   linear function that reduces phage burst size on susceptibles with respect to bacterial growth rate
    return args.phage_burstsize_linear*growthrate(nutrients)+args.phage_burstsize_const


def f(y,yd):
  # unpack arrays
  n,s,si,r,ri,p = y
  spDELAY,rpDELAY = yd
  # assume phage concentration is contant at boundaries
  # shifted arrays necessary for diffusion term of phage movement
  pfwd  = np.concatenate((p[1:],np.array([p[-1]])))
  pback = np.concatenate((np.array([p[0]]),p[:-1]))
  
  gr = growthrate(n)

  return np.array([	-gr*args.nutrients_per_cell*(s+r), \
			gr*s - args.phage_bacteria_absorptionconstant*s*p, \
			args.phage_bacteria_absorptionconstant *(s*p-spDELAY), \
			gr*r - args.phage_bacteria_absorptionconstant*args.bacteria_resitant_absorptionreduction*p*r, \
			args.phage_bacteria_absorptionconstant*args.bacteria_resitant_absorptionreduction*(r*p - rpDELAY) ,\
			burstsize(n)*args.phage_bacteria_absorptionconstant*(spDELAY+args.bacteria_resitant_absorptionreduction*rpDELAY)-args.phage_bacteria_absorptionconstant*(s+si+r+ri)*p + args.phage_diffusionconstant/args.dx**2*(pback - 2*p + pfwd)])

def RungeKutta4(y,yd):
  # 4th order Runge-Kutta integration scheme
  k1 = args.epsilon*f(y,yd[0])
  k2 = args.epsilon*f(y+k1/2.,(yd[0]+yd[1])/2.)
  k3 = args.epsilon*f(y+k2/2.,(yd[0]+yd[1])/2.)
  k4 = args.epsilon*f(y+k3,yd[1])
  ret = y + (k1+2*k2+2*k3+k4)/6.

  if args.cutatzero:
      ret[ret<0] = 0
  else:
    # restrict solutions for cells/phages to have more than one cell/phage in each bin
    (ret[1:])[ret[1:]<1] = 0
    # nutrient concentration should be larger than 0
    (ret[0:])[ret[0:]<0] = 0
  return ret

def main():
    parser = argparse.ArgumentParser()
    
    # different argument groups allow additional information to be printed to the screen
    # in the end everything is still parsed to the same "args" variable
    
    # define lattice
    parser_lattice = parser.add_argument_group(description="Lattice parameters")
    parser_lattice.add_argument("-s","--space",type=int,default=100)
    parser_lattice.add_argument("-d","--dx",type=float,default=0.2)
    
    # algorithm parameters
    parser_algorithm = parser.add_argument_group(description="Algorithm parameters")
    parser_algorithm.add_argument("-e","--epsilon",type=float,default=1e-4) # measured in hrs, thus default eps=1/200h=0.3min
    parser_algorithm.add_argument("-M","--maxsteps",type=int,default=300000)
    parser_algorithm.add_argument("-O","--outputsteps",type=int,default=500)
    parser_algorithm.add_argument("-o","--outputbins",type=int,default=1)
    parser_algorithm.add_argument("-Q","--quiet",default=False,action="store_true")
    parser_algorithm.add_argument("-Z","--cutatzero",default=False,action="store_true")

    # interactions between various trophic levels
    parser_interactions = parser.add_argument_group(description="Interactions between different trophic levels")
    parser_interactions.add_argument("-D","--phage_diffusionconstant",type=float,default=2.5e-2)
    parser_interactions.add_argument("-b","--phage_burstdelay",type=float,default=.5)
    parser_interactions.add_argument("-w","--phage_delaydistr",type=float,default=.05)
    parser_interactions.add_argument("-L","--phage_burstsize_linear",type=float,default=100)
    parser_interactions.add_argument("-l","--phage_burstsize_const",type=float,default=2)
    parser_interactions.add_argument("-A","--phage_bacteria_absorptionconstant",type=float,default=1e-6)
    parser_interactions.add_argument("-m","--bacteria_growth_Kmax",type=float,default=0.7204)
    parser_interactions.add_argument("-c","--bacteria_growth_Kc",type=float,default=0.0000257024)
    parser_interactions.add_argument("-n","--nutrients_per_cell",type=float,default=1e-10)          # mg/cell
    parser_interactions.add_argument("-r","--bacteria_resitant_absorptionreduction",type=float,default=1e-3)

    # startconcentrations 
    parser_startconc = parser.add_argument_group(description="Start concentrations")
    parser_startconc.add_argument("-S","--startconcentrationsusceptibles",type=float,default=10e4) # Indiv/mm
    parser_startconc.add_argument("-R","--startconcentrationresistant",   type=float,default=0e4)  # Indiv/mm
    parser_startconc.add_argument("-P","--startconcentrationphages",      type=float,default=100)  # Indiv in first bin
    parser_startconc.add_argument("-N","--startconcentrationnutrients",   type=float,default=1e-4) # mg/mm
    parser_startconc.add_argument("-p","--initialplaquesize",             type=float,default=.5)   # radius in mm


    global args
    args = parser.parse_args()
    
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
    # not concentrations, but rather numbers of cells/phages in a given bin
    y[1,:] = args.startconcentrationsusceptibles
    y[2,:] = 0
    y[3,:] = args.startconcentrationresistant
    y[4,:] = 0
    y[5,:int(args.initialplaquesize/args.dx)] = args.startconcentrationphages

    for o in range(args.maxsteps):
	# delay as convolution with 
	yd = np.array([[np.dot(delaydistr1,sdelay),np.dot(delaydistr1,rdelay)],[np.dot(delaydistr2,sdelay),np.dot(delaydistr2,rdelay)]])
	
	# integrate system of differential equations with Runge-Kutta method
	y = RungeKutta4( y ,yd)
 	
	# update stored solutions:
	# drop first item (0 as python index), append new solutions from current timestep to the end
	sdelay[0:delaysize-1,:] = sdelay[1:delaysize,:]
	sdelay[delaysize-1,:] = y[1]*y[5]
	rdelay[0:delaysize-1,:] = rdelay[1:delaysize,:]
	rdelay[delaysize-1,:] = y[3]*y[5]

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
  
  
  

  