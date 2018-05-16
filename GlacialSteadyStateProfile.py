#gprasicek 2017

import numpy as np
import matplotlib.pyplot as plt


def dhdxfromFlux(q,fd=7.26e-5,fs=3.27,us=1):
    # calculate ice thickness and ice sueface slope
    a = fd/fs*us
    c = us
    d = -q
    delta0 = -3*a*c
    delta1 = 27*a**2*d
    
    C = ((delta1 + (delta1**2-4*delta0**3)**(0.5))/2)**(1./3)
    
    H = -(1/(3*a))*(C+delta0/C)

    dhdx = (us/fs)**(1./n)*H**((1-n)/n)

    return [dhdx,H]

# set x domain and constants
length = 5e4
nx = 5e3
dx = length/(nx-1)
Exfrac = 0.5 # set relative horizontal position of equilibrium line
Ex = length*Exfrac
Expos = int(nx*Exfrac)

betaX = -5e-5
beta = 1e-3
K = 1e-4 
U = 1e-3
l = 1.
n = 3.
ELA = 0

ro = 910.0
g = 9.81 
yr = 31556926.0
Ad = 2.39e-24
As = 1.46e-19
fd = Ad*(ro*g)**n*yr
fs = As*(ro*g)**n*yr

# calculate steady state  
x = np.arange(1,nx-3)*dx
us = (U/K)**(1./l)
B = betaX*(x-length*Exfrac)
MBerr = 1e6        
        
while MBerr > 1:
    
    q = np.cumsum(B)*dx

    dhdx, H = dhdxfromFlux(q,fd,fs,us)
    
    h = np.cumsum(dhdx)*dx
    h = h.max()-h
    h = h+(ELA-h[Expos])
    
    Bideal = beta*(h-ELA)

    MBerr = np.sum(((B-Bideal)**2)**0.5)
    
    B = Bideal
    

Zb = h-H

# Plot

plt.xlabel("Distance [m]")
plt.ylabel("Elevation [m]")
        
plt.plot(x, h, 'k-', label='Ice surfae')
plt.plot(x, Zb, 'k--', label='Bedrock')

plt.legend(loc=1, fontsize=9, frameon=False)

