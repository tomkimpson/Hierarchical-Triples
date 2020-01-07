from __future__ import division
import numpy as np
from setup import *
from derivatives import *
import sys
from EvalF import Fn
import time


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
    print ('Numerical orbital evolution has completed with fs = ', fs, ' Hz and Tobs = ', Tint, ' years')
    return output,const

def analytical_orbital_evolution(fit_data,Tint,fs,const):
    
    #Extract the relevant constants
    #K, J2, C,A,eta,I,Lambda
    K = const[0]
    J2 = const[1]
    C = const[2]
    I = const[5]
    lam = const[6]

    
    #time
    T_seconds = Tint*365*24*3600
    t = np.arange(0,T_seconds,1/fs)
    #t = fit_data[:,0]
    

    #Get the eccentricity behaviour
    tstart = time.time()
    A,B,omega,offset = extract(fit_data,1)
    tend = time.time()
    print ('The eccentricity fit completed in', tend-tstart,'seconds')
    e_approx = A*np.sin(omega*t) +B*np.cos(omega*t)  + offset 

    

    

    
    #Semi major axis
    F0 = Fn(A,B,omega,offset,0)
    Fbar = Fn(A,B,omega,offset,t)
    
    normalisation = fit_data[0,3]**4 / 4 - C*F0
    
    a_approx = (4*C*Fbar + 4*normalisation)**(1/4)

    
    #Precession of periastron
    #---First get the analytical derivative of gamma
    t1 = fit_data[:,0]
    e1 = fit_data[:,1]
    g1 = fit_data[:,2]
    a1 = fit_data[:,3] 
    J1 = fit_data[:,4] 
    u = 1 - e1**2
        
    NSP = 2*K*a1**2 * u * (2/J1 + np.cos(I)/J2) + lam * a1**(-2.5) * u**(-1)
    SP = 10*K*a1**2*(e1**2*np.cos(I)/J2 - (u-np.cos(I)**2)/J1)
    
    gderiv = NSP + SP*np.sin(g1)**2
    
    print(len(gderiv) ,len(t), len(fit_data))
    
    #---Now fit the derivative with a parametric function
    B_A, B_omega, B_offset, B_D,gamma_f, gamma_g,H = extract2(t1,gderiv)
    #We can directly integrate this function
    integration_constant = g1[0] - integral(B_A, B_omega, B_offset, B_D,gamma_f, gamma_g,H,0) 
    g_approx = integral(B_A, B_omega, B_offset, B_D,gamma_f, gamma_g,H,t) + integration_constant
    

    

    print (g_approx)
    
    
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
    print ('extract ecc')
    
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


def AdvancedTrig(t,B_A, B_omega, B_offset, B_D,gamma_f, gamma_g,H):
    B = B_A*np.sin(B_omega*t + B_offset) +B_D
    s2g = np.sin(gamma_f*t + gamma_g)**2
    return H + B*s2g
    
def extract2(t,f):
    print ('Extract for double trig func')
    func = AdvancedTrig
    
    H= 3.631335880308346e-05
    B_A = 2.50885879e-05
    B_omega = 4.32315349e-05
    B_offset = -1.76843846e-01 
    B_D = -2.17565784e-05
    gamma_f = 2.17209725e-05
    gamma_g = 6.85614205e-01
    
    
    
    #H= 1e-5
   # B_A = 1e-5
   # B_omega = 1e-5
   # B_offset = 1e-1
   # B_D = 1e-5
   # gamma_f = 1e-5
   # gamma_g = 1e-1
    

    p0 = (B_A, B_omega, B_offset, B_D,gamma_f, gamma_g,H)   
    popt, pcov = curve_fit(func, t,f,p0=p0)

    #return func(t, *popt)
    return popt


def integral(A,omega,C,D,f,g,H,t):
    
    part1 = -A*np.cos(C+t*(omega-2*f)-2*g)/(4*(2*f-omega))
    
    part2 =  A*np.cos(C+t*(omega+2*f)+2*g)/(4*(2*f+omega))
    
    part3 = 0.5*A*(np.sin(C)*np.sin(t*omega)/omega - np.cos(C)*np.cos(omega*t)/omega)
    
    part4 = D*(f*t + g)/(2*f)
    
    part5 = -D*np.sin(2*(f*t + g))/(4*f)
    
    part6 = H*t
    
   
    return part1 + part2 +part3 + part4 + part5 + part6

# --- TESTING---