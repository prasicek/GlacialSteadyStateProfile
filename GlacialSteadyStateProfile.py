#gprasicek 2017

import numpy as np
import matplotlib.pyplot as plt


def ISSfromFlux(q,fd=7.26e-5,fs=3.27,us=1):
    # calculate ice thickness and ice sueface slope
    a = fd/fs*us
    c = us
    d = -q
    delta0 = -3*a*c
    delta1 = 27*a**2*d
    
    C = ((delta1 + (delta1**2-4*delta0**3)**(0.5))/2)**(1./3)
    
    H = -(1/(3*a))*(C+delta0/C)

    S = (us/fs)**(1./3)*H**(-2./3)

    return [S,H]

# set x domain and constants
length = 5e4
lengthfrac = 0.5 # set relative horizontal position of equilibrium line

betaX = -5e-5
beta = 1e-3
K = 1e-4 
U = 1e-3
l = 1
ELA = 0

fs = 3.27
fd = 5.388e-5

# calculate steady state  
x = np.arange(length)
us = (U/K)**(1./l)
acc = betaX*(x-length*lengthfrac)
MBerr = 1e6        
        
while MBerr > 1:
    
    q = np.cumsum(acc) 

    S, H = ISSfromFlux(q,fd,fs,us)
    S[np.isnan(S)] = 0
    
    Zs = np.cumsum(S)
    Zs = Zs.max()-Zs
    ELApos = np.abs(acc).argmin()
    Zs = Zs+(ELA-Zs[ELApos])
    
    accideal = beta*(Zs-ELA)

    MBerr = np.sum(((acc-accideal)**2)**0.5)
    
    acc = accideal
    
Zs = Zs[H>0]
acc = acc[H>0]
q = q[H>0]
S = S[H>0]   
x = x[H>0]
H = H[H>0]

Zb = Zs-H

# Plot

plt.xlabel("Distance [m]")
plt.ylabel("Elevation [m]")
        
plt.plot(x, Zs, 'k-', label='Ice surfae')
plt.plot(x, Zb, 'k--', label='Bedrock')

plt.legend(loc=1, fontsize=9, frameon=False)

