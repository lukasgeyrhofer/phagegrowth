#!/usr/bin/env python

import numpy as np
import argparse
import sys,math


def RungeKutta4(func,y,epsilon):
  # 4th order Runge-Kutta integration scheme
  k1 = epsilon*func(y)
  k2 = epsilon*func(y+k1/2.)
  k3 = epsilon*func(y+k2/2.)
  k4 = epsilon*func(y+k3)
  ret = y + (k1+2*k2+2*k3+k4)/6.
  ret[ret<0] = 0.
  return ret


def phagedynamics(y):
    if y[3] <= 0:
        growthrate = 0
    else:
        growthrate = param['growthrate']
        
    return np.array([ (growthrate - y[1]*y[2]*param['absorption']),
                     -param['absorption']*y[1]*(1-y[1])*y[2]/y[0],
                     (param['burstsize']*y[1]-1)*param['absorption']*y[2],
                     -growthrate*   param['invyield']*y[0]])

def output(time,concentrations,widthtime = 4):
    print "{value:.{width}f}".format(value=time,width = widthtime),
    for c in concentrations:
        print "{:.5e}".format(c),
    print

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-B","--initial_bacteria",type=float,default=1e5)
    parser.add_argument("-S","--initial_susceptible_fraction",type=float,default=.5)
    parser.add_argument("-P","--initial_phage",type=float,default=1e1)
    parser.add_argument("-N","--initial_nutrients",type=float,default=1e-2)
    
    parser.add_argument("-a","--param_growthrate",type=float,default=.7)
    parser.add_argument("-b","--param_burstsize",type=float,default=100)
    parser.add_argument("-n","--param_absorption",type=float,default=0.007343)
    parser.add_argument("-y","--param_nutrientpercell",type=float,default=1e-10)
    
    parser.add_argument("-e","--algorithm_epsilon",type=float,default=1e-3)
    parser.add_argument("-T","--algorithm_maxtime",type=float,default=10)
    parser.add_argument("-O","--algorithm_outputstep",type=int,default=100)
    
    args = parser.parse_args()
    
    global param
    param = {'growthrate' : args.param_growthrate, 'burstsize': args.param_burstsize, 'absorption': args.param_absorption, 'invyield': args.param_nutrientpercell}
    
    y = np.array([args.initial_bacteria,args.initial_susceptible_fraction, args.initial_phage, args.initial_nutrients])
    
    i = 0
    maxsteps = args.algorithm_maxtime / args.algorithm_epsilon
    outputwidth = int(np.rint(-np.log10(args.algorithm_epsilon * args.algorithm_outputstep)))
    print outputwidth
    output(0,y,outputwidth)
    while i < maxsteps:
        y = RungeKutta4(phagedynamics,y,args.algorithm_epsilon)
        i += 1

        if i % args.algorithm_outputstep == 0:
            output(i*args.algorithm_epsilon,y,outputwidth)


if __name__ == "__main__":
    main()
    
    

