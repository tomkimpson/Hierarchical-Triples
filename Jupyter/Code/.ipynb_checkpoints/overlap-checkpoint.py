from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


def overlap(data1,data2):
    
    t = data1[:,0]
    hplus = data1[:,1]
    hcross = data1[:,2]
    
    
    dt = t[1] - t[0]
    fs = 1/dt
   

    #Get the frequencies
    f = np.fft.rfftfreq(t.size,dt)

    

    #Take the FT
    hp = dt*np.fft.rfft(hplus)
    hc = dt*np.fft.rfft(hcross)
    
    hN =np.sqrt(abs(hp)**2 + abs(hc)**2) 
    
    
    #Get rid of zeroth frequencies - WHY?
    hp = hp[1:] # get rid of zeroth frequency
    hc = hc[1:]
    f = f[1:]
    
    #Noise
    print (min(f),max(f),len(f))
    #f = np.linspace(1e-6,1,10000)
    S = noise(f)
    
 

    #output
    overlap = 0.50
    hsig = np.sqrt(hp**2 + hc**2)
    
    return f,hsig, np.sqrt(S)

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
def noise(f):
    #Calculate the LISA noise curve
    Larm = 2.5e9
    Clight = 3e8
    fstar = Clight/(2*np.pi*Larm)
    NC = 2

    alpha = 0.133
    beta = 243.
    kappa = 482.
    gamma = 917.
    f_knee = 2.58e-3

    A = 1.8e-44/NC
    Sc = 1. + np.tanh(gamma*(f_knee-f))
    Sc *=np.exp(-f**alpha + beta*f*np.sin(kappa*f))

    Sc *= A*f**(-7./3.)


    #LISA response function

    RFILE = np.loadtxt('ResponseFunction.txt')
    Rx = RFILE[:,0] * fstar
    Ry = RFILE[:,1] * NC

    newR = np.interp(f,Rx,Ry)
    R = newR


    #Power Spectral Density
    P_oms = (1.5e-11)**2 * (1. + (2.0e-3/f)**4)
    P_acc = (3.0e-15)**2 * (1.+(0.4e-3/f)**2)*(1. + (f/(8.0e-3))**4)
    Pn = (P_oms + 2.*(1. + np.cos(f/fstar)**2)*P_acc/(2*np.pi*f)**4)/Larm**2

    #Total noise
    S = Pn/R + Sc
    

    
    return S