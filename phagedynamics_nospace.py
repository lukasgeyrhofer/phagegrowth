#!/usr/bin/env python
# -*- coding: utf-8 -*-

# dynamics of phages without spatial component
# numerical solutions based on Runge-Kutta integration (of 4th order)
#
# ========================================================
#   set of equations
#   bacterial populations are Bx, infected bacteria are Ix, phage is P, nutrients are S
# ========================================================
#
#    d/dt B0 = a B0 - e B0 P
#    d/dt B1 = a B1
#    d/dt I0 = e B0 P - 1/t (f * I0)
#    d/dt I1 = 0
#    d/dt P  = b/t (f * I0) - eta * (B0 + B1 + I0 + I1) * P
#    d/dt S  = -a/Y (B0 + B1 + I0 + I1)
#
# =======================================================
#   convolution with Gaussian burst delay kernel is denoted as "f *"
#   growth rate a depends on nutrient concentration via Monod growth, a = amax S/(Km + S)
#   varying phage growth is absorbed as linear interpolation into burstsize b, b = (Bmax-Bmin) a/amax + Bmin
# =======================================================

import numpy as np
import argparse
import sys,math


def RungeKutta4(func,y,yd,epsilon):
  # 4th order Runge-Kutta integration scheme
  k1 = epsilon*func(y,yd[0])
  k2 = epsilon*func(y+k1/2.,np.mean(yd,axis=1))
  k3 = epsilon*func(y+k2/2.,np.mean(yd,axis=1))
  k4 = epsilon*func(y+k3,yd[1])
  ret = y + (k1+2*k2+2*k3+k4)/6.
  ret[ret<0] = 0.
  return ret



def phagedynamics(y,yd):
    # monod kinetics for bacterial growth
    growthrate = param['growthrate'] * y[5]/(param['growthkm'] + y[5])
    # define burstsize as compound parameter for phage growth, depends linearly on bacterial growthrate and interpolates between minimal and maximal burstsizes
    burstsize = growthrate/param['growthrate'] * ( param['phagemaxburst'] - param['phageminburst']) + param['phageminburst']
    
    return np.array([growthrate * y[0] - param['absorption']*y[0]*y[4],                             # 0: susceptible bacteria
                     growthrate * y[1],                                                             # 1: resistant bacteria
                     param['absorption'] * y[0] * y[4] - yd[0],                                     # 2: susceptible bacteria infected
                     0,                                                                             # 3: resistant bacteria infected
                     burstsize * yd[0] - param['absorption'] * (y[0] + y[1] + y[2] + y[3]) * y[4],  # 4: phage
                     -growthrate*param['invyield']*(y[0] + y[1] + y[2] + y[3]) ])                   # 5: nutrients
    
        
        
def output(time,concentrations,widthtime = 4):
    print "{value:.{width}f}".format(value=time,width = widthtime),
    for c in concentrations:
        print "{:12.5e}".format(c),
    print



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-B","--initial_bacteria",type=float,default=1e5)
    parser.add_argument("-S","--initial_susceptible_fraction",type=float,default=.5)
    parser.add_argument("-P","--initial_phage",type=float,default=1e1)
    parser.add_argument("-N","--initial_nutrients",type=float,default=1e-4)
    
    parser.add_argument("-a","--param_growthrate",type=float,default=.7)
    parser.add_argument("-b","--param_burstsize",type=float,default=100)
    parser.add_argument("-n","--param_absorption",type=float,default=1e-4)
    parser.add_argument("-y","--param_nutrientpercell",type=float,default=1e-10)
    parser.add_argument("-k","--param_growth_km",type=float,default = 0.0000257024)
    parser.add_argument("-g","--param_minburst",type=float,default=2)
    parser.add_argument("-G","--param_maxburst",type=float,default=100)
    
    parser.add_argument("-t","--timedelay_mean",type=float,default=0.5)     # latent period is half hour
    parser.add_argument("-s","--timedelay_stddev",type=float,default=0.08)  # with a standard deviation of roughly 5min
    
    parser.add_argument("-e","--algorithm_epsilon",type=float,default=1e-3)
    parser.add_argument("-T","--algorithm_maxtime",type=float,default=20)   # in hrs
    parser.add_argument("-O","--algorithm_outputstep",type=int,default=100)
    
    args = parser.parse_args()

    global param
    param = {'growthrate' : args.param_growthrate, 'burstsize': args.param_burstsize, 'absorption': args.param_absorption, 'invyield': args.param_nutrientpercell, 'growthkm': args.param_growth_km, 'phageminburst': args.param_minburst, 'phagemaxburst':args.param_maxburst, 'timedelay_mean' : args.timedelay_mean, 'timedelay_stddev' : args.timedelay_stddev}
    outputwidth = int(np.ceil(-np.log10(args.algorithm_epsilon * args.algorithm_outputstep))) + 1
    maxsteps = int(args.algorithm_maxtime / args.algorithm_epsilon)
    

    # store solutions to implement time-delay in phage burst
    delaysize = int((args.timedelay_mean+3*args.timedelay_stddev)/args.algorithm_epsilon)
    delaytime = np.arange(delaysize)*args.algorithm_epsilon-3*args.timedelay_stddev
    # two distributions shifted by a timestep for use in Runge-Kutta algorithm
    delaydistr1 = np.exp(-0.5*(delaytime/args.timedelay_stddev)**2)/np.sqrt(2*math.pi*args.timedelay_stddev**2)*args.algorithm_epsilon
    delaydistr2 = np.exp(-0.5*((delaytime+args.algorithm_epsilon)/args.timedelay_stddev)**2)/np.sqrt(2*math.pi*args.timedelay_stddev**2)*args.algorithm_epsilon
    
    # generate initial conditions
    y = np.array([args.initial_bacteria * args.initial_susceptible_fraction, args.initial_bacteria* (1. - args.initial_susceptible_fraction), 0., 0., args.initial_phage, args.initial_nutrients])
    infecteddelay = np.zeros((delaysize,2))
    yd = np.array([np.dot(delaydistr1,infecteddelay),np.dot(delaydistr2,infecteddelay)])
    
    output(0,y,outputwidth)

    for i in range(1,maxsteps+1):

        infecteddelay[0:delaysize-1,:] = infecteddelay[1:delaysize,:]
        infecteddelay[delaysize-1,0] = y[2]
        infecteddelay[delaysize-1,1] = y[3]
        yd = np.array([np.dot(delaydistr1,infecteddelay),np.dot(delaydistr2,infecteddelay)])

        y = RungeKutta4(phagedynamics,y,yd,args.algorithm_epsilon)
        
        if i % args.algorithm_outputstep == 0:
            output(i*args.algorithm_epsilon,y,outputwidth)


if __name__ == "__main__":
    main()
    
    

