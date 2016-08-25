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
    
    bact = y[0] + y[1] + y[2] + y[3]
    
    if y[-1] <= 0: # nutrients are depleted
        growthrate  = 0
        burstsize   = param['burstsize_depletion']
        latencytime = param['latencytime_depletion']
    else:
        growthrate  = param['growthrate']
        burstsize   = param['burstsize']
        latencytime = param['latencytime']
    
    # new model for recovery
    # assume linear relation F(MOI) from our Fig4 'Burst probability of resistants', and set
    # b/(r+b) = F(MOI)
    # where b and r are burst and recovery rates
    #recovery = 0
    #if bact > 0:
        #recovery = 1/latencytime * (param['resistant_maxMOI'] - y[4]/bact)/(param['resistancereductionrate'] + y[4]/bact)
    #if recovery < 0:
        #recovery = 0
    # dont have enough experimental estimates to properly infer parameters
    # fall back to simpler assumption that MOI does not influence the recovery
    recovery = param['resistancereductionrate']
    
    
    return np.array([   growthrate * y[0] - param['absorption'] * y[0] * y[4],                      # 0  susceptible bacteria
                        growthrate * y[1] - param['absorption'] * y[1] * y[4] + recovery * y[3],    # 1  resistant bacteria
                        param['absorption'] * y[0] * y[4] - y[2]/latencytime,                       # 2  infected susceptible bacteria
                        param['absorption'] * y[1] * y[4] - y[3]/latencytime - recovery*y[3],       # 3  infected resistant bacteria
                        burstsize/latencytime * (y[2] + [3]) - param['absorption']*y[4]*bact,       # 4  phages
                        -growthrate*param['nutrientspercell']*(y[0] + y[1])   ])                    # 5  nutrients

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
    parser.add_argument("-N","--initial_nutrients",type=float,default=1.) # in dilutions of the original LB medium
    
    parser.add_argument("-a","--param_growthrate",type=float,default=.63)
    parser.add_argument("-b","--param_burstsize",type=float,default=90.)
    parser.add_argument("-l","--param_latencytime",type=float,default=.5)
    parser.add_argument("-d","--param_burstsize_depletion",type=float,default=3)
    parser.add_argument("-L","--param_latencytime_depletion",type=float,default=1.8)
    parser.add_argument("-n","--param_absorption",type=float,default=1e-7)
    parser.add_argument("-y","--param_nutrientspercell",type=float,default=2e-10) # also in dilutions of original LB medium
    parser.add_argument("-r","--param_resitant_reduction_rate",type=float,default=2.)
    #parser.add_argument("-R","--param_resitant_maxMOI",type=float,default=1e6)
    
    parser.add_argument("-e","--algorithm_epsilon",type=float,default=1e-3)
    parser.add_argument("-T","--algorithm_maxtime",type=float,default=48)
    parser.add_argument("-O","--algorithm_outputstep",type=int,default=100)
    
    args = parser.parse_args()
    
    global param
    param = {'growthrate' : args.param_growthrate, 'burstsize': args.param_burstsize, 'burstsize_depletion': args.param_burstsize_depletion, 'latencytime':args.param_latencytime,'latencytime_depletion':args.param_latencytime_depletion, 'absorption': args.param_absorption, 'nutrientspercell': args.param_nutrientspercell,'resistancereductionrate': args.param_resitant_reduction_rate}
    
    y = np.array([args.initial_bacteria*args.initial_susceptible_fraction,args.initial_bacteria*(1. - args.initial_susceptible_fraction),0.,0., args.initial_phage, args.initial_nutrients])
    
    i = 0
    maxsteps = args.algorithm_maxtime / args.algorithm_epsilon
    outputwidth = int(np.rint(-np.log10(args.algorithm_epsilon * args.algorithm_outputstep)))
    output(0,y,outputwidth)
    while i < maxsteps:
        y = RungeKutta4(phagedynamics,y,args.algorithm_epsilon)
        i += 1

        if i % args.algorithm_outputstep == 0:
            output(i*args.algorithm_epsilon,y,outputwidth)


if __name__ == "__main__":
    main()
    
    

