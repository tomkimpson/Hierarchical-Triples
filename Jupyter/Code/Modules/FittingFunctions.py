from scipy.optimize import curve_fit 
import numpy as np


def FitEcc(x,y):
    
    #Define the function fit
    func = ansatz
    
    #Define the initial guess of parameters
    #OSC
    A,B,omega = 0.05,0.05,1e-4

    #TAPER
    k_taper = 2e-7
    
    #DECAY
    k_decay = 1e-6
    G = 1/(np.exp(k_decay*x[-1] -1))
    F = 1 + G

    decay = F-G*np.exp(k_decay*x)
    
    
    
    sigma = np.ones(len(x))
    #sigma[[-1]] = 1e-1 #weight the end point heavily
    

    
    #SHIFT
    delta = 0.9


    #Put into vector
    p0 = (A,B,omega,
          k_taper,
         F,G,k_decay,
         delta)
    

    #Do the fit
    popt, pcov = curve_fit(func, x,y,p0=p0,sigma=sigma)
    
    
    output = np.zeros((len(x),2))
    output[:,0] = x
    output[:,1] = func(x, *popt)
    
    print (popt)
    
    return output
    
    

    
    
    
def ansatz(x,A,B,omega,
           k_taper,
           F,G,k_decay,
           delta
          ):
    
    #oscillatory part
    osc = A*np.sin(omega*x) + B*np.cos(omega*x) 
    
    #Taper function
    taper = np.exp(-k_taper*x)
    
    #decay function
    decay = F-G*np.exp(k_decay*x)


    
    return  taper*osc+decay*delta
    
    
    
