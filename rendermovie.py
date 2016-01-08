#!/usr/bin/env python

import numpy as np
import argparse
import sys,math
import Gnuplot

parser = argparse.ArgumentParser()
parser.add_argument("-c","--conffile",required=True)
parser.add_argument("-N","--normalize",action='store_true',default=False)
parser.add_argument("-m","--maxlength",type=float,default=None)
args = parser.parse_args()

try:
    data = np.genfromtxt(args.conffile)
except:
    print >> sys.stderr,"could not open file"
    exit(1)

tall  = data[:,0]
xall  = data[:,1]
nall  = data[:,2]
sall  = data[:,3]
siall = data[:,4]
rall  = data[:,5]
riall = data[:,6]
pall  = data[:,7]

tlist = np.unique(tall)

dt = tall[1] - tlist[0]

x = xall[tall == tall[1]]
if args.maxlength == None:
    xmax = x[-1]
else:
    xmax = args.maxlength


m = np.ones(6)
if args.normalize:
    m[0] = max(nall)
    m[1] = max(sall)
    m[2] = max(siall)
    m[3] = max(rall)
    m[4] = max(riall)
    m[5] = max(pall)

frame = 0


print m
#exit(1)

g = Gnuplot.Gnuplot(persist=0)


for t in tlist:
    print frame, t
    
    n = nall[tall==t]/m[0]
    s = sall[tall==t]/m[1]
    si = siall[tall==t]/m[2]
    r = rall[tall==t]/m[3]
    ri = riall[tall==t]/m[4]
    p = pall[tall==t]/m[5]
    
    #print x,n,s

    


    pn  = Gnuplot.Data(x,n ,with_="lines lw 2 lc rgb \"#f57900\"", title="N(t,x)")
    ps  = Gnuplot.Data(x,s ,with_="lines lw 2 lc rgb \"#204a87\"", title="S(t,x)")
    psi = Gnuplot.Data(x,si,with_="lines lw 2 lc rgb \"#729fcf\"", title="SI(t,x)")
    pr  = Gnuplot.Data(x,r ,with_="lines lw 2 lc rgb \"#4e9a06\"", title="R(t,x)")
    pri = Gnuplot.Data(x,ri,with_="lines lw 2 lc rgb \"#8ae234\"", title="RI(t,x)")
    pp  = Gnuplot.Data(x,p ,with_="lines lw 4 lc rgb \"#ef2929\"", title="P(t,x)")

    g.reset()
    
    g("set terminal pngcairo enhanced size 1000,400")
    g("set output \"frame%04d.png\""%frame)
    g("unset key")
    g("set border 15 lw 2")
    g("set tics front")
    g("set key bottom right spacing 2")
    frame += 1
    
    g("set xra [0:%lf]"%(xmax))
    #g("set xtics 500")
    
    if args.normalize:
        g("set yra [0:1.1]")
        g("set ytics (.5,\"max\"1)")
    else:
        g("set yra [0:%lf]"%1.1*max(m))
    
    g("set label 1 \"t = %.1lf\" at screen .25,.9 center front"%t)
    
    g.plot(pn,ps,psi,pr,pri,pp)
    
    del pn
    del ps
    del pr
    del pp
    del psi
    del pri

