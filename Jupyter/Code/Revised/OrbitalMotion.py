from __future__ import division
import numpy as np
from setup import *
from derivatives import *
import sys
from EvalF import Fn



def univeral_constants():
    c = 3e8
    G = 6.67e-11
    Msolar = 2e30
    AU = 1.49e11 #meters
    pc = 3.086e16 #parsec in m
    
    return c,G,Msolar,AU, pc


#----PUBLIC ----
def numerical_orbital_evolution(m0,m1,m2,a1,e1,g1,J1,e2,a2,I,fs,Tint):
    
    
    #Get some constants of the motion
    const = constants(m0,m1,m2,e2,a2,I)
    
    yn = np.array((e1,g1,a1,J1))  
    
    #Call RK solver
    output = RungeKutta(yn,const,fs,Tint)
    
    #output
    #output[:,2] = np.sin(2*output[:,2]) #Convert from gamma to sin 2 gamma
    print ('Numerical orbital evolution has completed with fs = ', fs, ' Hz and Tobs = ', Tint, ' years')
    return output

def analytical_orbital_evolution(fit_data,Tint,fs):
    
    #Get some constants of the motion
    #const = constants(m0,m1,m2,e2,a2,I)
    C = -6.752505469155556e+23
    
    #time
    T_seconds = Tint*365*24*3600
    t = np.arange(0,T_seconds,1/fs)
    #t = fit_data[:,0]
    

    #Get the eccentricity behaviour
    A,B,omega,offset = extract(fit_data,1)
    e_approx = A*np.sin(omega*t) +B*np.cos(omega*t)  + offset #+ delta*t
    
    #Semi major axis
    

    F0 = Fn(A,B,omega,offset,0)
    Fbar = Fn(A,B,omega,offset,t)
    
    normalisation = fit_data[0,3]**4 / 4 - C*F0
    
    a_approx = (4*C*Fbar + 4*normalisation)**(1/4)
    #a_approx = np.ones(len(t)) * fit_data[0,3] #constant GW
    
    #Precession of periastron
    g_approx = extract2(fit_data,2,t) 
    
    
    
    #------- failed attempt at foueri decom - too computationally challenging
    
    
    #xx = fit_data[:,0]
   # yy = fit_data[:,2]
   # period = 2*np.pi/omega
    

   # freqs = np.fft.fftfreq(len(xx), 1/fs)
   # Fk = np.ones(len(xx))
    
   # for k in range(len(Fk)):
   #     for j in range(len(yy)):
            
   #         ar = yy[j]*np.exp(-1j * 2*np.pi * freqs[k]*yy[j])/fs
   #         Fk[k] = np.sum(ar)
        
        
   # fx = np.ones(len(xx))
   # for k in range(len(xx)):
    #    for j in range(len(Fk)):
    #        ar = Fk[j] * np.exp(1j * 2*np.pi *freqs[j]*xx[k])*fs
    #        fx[k] = np.sum(ar)
    
    
   
  
    
  

    
    
    #output
    out = np.ones((len(t),5))
    out[:,0] = t
    out[:,1] = e_approx
    out[:,2] = g_approx
    out[:,3] = a_approx
    #out[:,4] = fit_data[:,4]



    
    #output
    print ('Analytical orbital evolution has completed with fs = ', fs, ' Hz and Tobs = ', Tint, ' years')
    return out

    
    
    
#----PRIVATE---



def RungeKutta(yn,const,fs,Tint):
    
    #Setup timing precision
    h = 1/fs
    Ts = Tint*365*24*3600
    trange = np.arange(0,Ts,h)
    tfinal = trange[-1]    
    
    t = 0
    nsteps = int(tfinal*fs) + 1

    
    #Define output array
    out = np.zeros((nsteps,5)) #t,e,gamma,a, J1
    counter = 0
    out[counter,0] = t
    out[counter,1] = yn[0]
    out[counter,2] = yn[1] 
    out[counter,3] = yn[2] 
    out[counter,4] = yn[3] 
    counter = counter + 1
    
    while t < tfinal:

        
        k1 = h * derivs(yn,const)
        k2 = h * derivs(yn+k1/2,const)
        k3 = h * derivs(yn+k2/2,const)
        k4 = h * derivs(yn+k3,const)

        ynew = yn + (k1 + 2*k2 + 2*k3 + k4)/6
        yn = ynew
        
   
    
        t = t + h
        
        if counter <= nsteps-1:
            out[counter,0] = t
            out[counter,1] = yn[0]
            out[counter,2] = yn[1]  
            out[counter,3] = yn[2] 
            out[counter,4] = yn[3] 
        
        
            counter = counter + 1    

    return out
    
    
    
def derivs(y,const):
    

    e_dot = edot(y,const)
    gamma_dot = gdot(y,const)
    a_dot = adot(y,const)
    J_dot = Jdot(y,const)


    return np.array((e_dot, gamma_dot,a_dot,J_dot))






def constants(m0,m1,m2,e2,a2,I):
    c,G,Msolar,AU, pc = univeral_constants()
    
    
    M = m0+m1
    mu = m0*m1/M
    
   
    
    
    K = 3*G*m0*m1*m2/(8*M) *a2**(-3) * (1-e2**2)**(-3/2)
    C = -64/5 * G**3*mu*M**2/c**5
    A = 5*K*np.sin(I)**2
    eta = -32/5 * G**(7/2) * mu**2 * M**(5/2) /c**5
    J2 = angular_momentum(M,m2,e2,a2)
    Lambda = 3* (G*M)**(3/2)/c**2
        
        
    print ('K = ', K)
    print ('C = ', C)
    print ('lambda = ', Lambda,G, c, m0,m1)
    print ('I = ', I)
    print ('J2 = ', J2)

        
    return np.array((K, J2, C,A,eta,I,Lambda)) 

    

    
    
from scipy.optimize import curve_fit   
def doube_trig_function(t,A,B,omega,offset):
    return A*np.sin(omega*t) + B*np.cos(omega*t)  + offset
    
def extract(data,index):
    t = data[:,0]
    f = data[:,index] 
    
    func = doube_trig_function
    offset = (f.max() + f.min()) / 2
    y_shifted = f - offset
    delta = 1e-6
    p0 = (
         (f.max() - f.min()) / 2,
         (f.max() - f.min()) / 2,
         np.pi * np.sum(y_shifted[:-1] * y_shifted[1:] < 0) / (t.max() - t.min()),
        offset
         )
        
    popt, pcov = curve_fit(func, t,f,p0=p0)
    return popt
    
    
 





    
    
# --- TESTING---  


def simple_trig(t,omega,offset):
    return omega*t + offset
    
def extract2(data,index,t1):
    t = data[:,0]
    f = data[:,index] 
    
    func = simple_trig
    offset = (f.max() + f.min()) / 2
    y_shifted = f - offset
    p0 = (
         np.pi * np.sum(y_shifted[:-1] * y_shifted[1:] < 0) / (t.max() - t.min()),
        offset
         )
        
    popt, pcov = curve_fit(func, t,f,p0=p0)
    print ('popt for gamma:', popt)
    return func(t1, *popt)







# --- TESTING---