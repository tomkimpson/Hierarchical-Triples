from __future__ import division
import numpy as np


def univeral_constants():
    c = 3e8
    G = 6.67e-11
    Msolar = 2e30
    AU = 1.49e11 #meters
    pc = 3.086e16 #parsec in m
    
    return c,G,Msolar,AU, pc




def status(fs, Tint):
    year = 365*24*3600
    T = Tint*year
    
    dt = 1/fs
    
    MaxF =fs/2
    MinF = 1 / (T)
    nsteps = T/dt
    
    print ('Sampling frequency = ', fs,  ' Hz for an observation period of ', Tint,' years')
    
    print ('Total number of integration steps is ', nsteps)
    print ('Frequency range is: ',  MinF,' - ', MaxF, ' Hz')

    
    
    
#Some basic Keplerian defenitions of third law, angular momentum etc.

def semi_major_axis(M,f):
    c,G,Msolar,AU, pc = univeral_constants()
    
    a = (G*M / (4*np.pi**2*f**2)) ** (1/3)
    
    return a
    
    
    
def angular_momentum(m0,m1,e,a):
    c,G,Msolar,AU, pc = univeral_constants()

    J = m0*m1*np.sqrt(G*a*(1-e**2)/(m0+m1))
    
    return J