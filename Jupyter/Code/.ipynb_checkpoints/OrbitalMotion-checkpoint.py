from __future__ import division
import numpy as np
from OrbitalDerivatives import *
c = 3e8
G = 6.67e-11
Msolar = 2e30
AU = 1.49e11 #meters



def get_orbital_evolution_play(m0,m1,f1,e1,beta,m2,e2,I,gamma,Tint,fs):
    
    print ('Getting the orbit')
    
    #Calculate some constants and useful values
    K,J1,J2,mu1,M,a1 = setup(m0,m1,m2,f1,e1,e2,beta)


    #Define time variable
    T_seconds = Tint*365*24*3600
    t = np.arange(0,T_seconds,1/fs)
    
    y = np.array((e1,gamma,a1))
    constants = np.array((K,J1,J2,I,mu1,M))
    
    
     #Get approx a

    da = adot(y,constants)
    approx_a = da*t + a1
    
    #MODS
    Cprime = -64*G**3 * mu1 * M**2 / (5 * c**5)
    ed = edot(y,constants)
    print ('edot =',ed)

    

    #Get approx gamma
    dg = gdot(y,constants)
    approx_gamma = t*dg + gamma #+ 0.01*np.sin(dg*t)
    
    
    print ('dg = ', dg)
    print ('da = ', da)
    
    
    #Approx e
    ge = e1*(1-e1**2)**(-5/2) * (1 + 121*e1**2/304)
    he = ge/(e1*(1-e1**2))
    psi = get_psi(approx_a,approx_gamma,da,dg,a1,gamma,t)
    AA = 5*K*(1-np.cos(I)**2)/J1
    Cprime = -64*G**3 * mu1 * M**2 / (5 * c**5)
    
    print ('Crpime = ',Cprime)

    CC = np.log(e1) - 0.5*np.log(1-e1**2)
    
    alpha_KL = AA*psi + CC
    alpha_GW = -19/36 * Cprime*he/da * (approx_a**(-3) - a1**(-3))

    alpha = alpha_KL + alpha_GW

    
    approx_e = np.exp(alpha)/np.sqrt(1+np.exp(2*alpha))

    
    #output
    out = np.zeros((len(t),4))
    out[:,0] = t
    out[:,1] = approx_e
    out[:,2] = approx_gamma
    out[:,3] = approx_a
    
    
    #out[:,0] = t
    #out[:,1] = e1
    #out[:,2] = gamma
    #out[:,3] = a1
    
    
    return out





def get_orbital_evolution_numerical(m0,m1,f1,e1,beta,m2,e2,I,gamma,Tint,fs):
     
    #Calculate some constants and useful values
    K,J1,J2,mu1,M,a1 = setup(m0,m1,m2,f1,e1,e2,beta)
    
    #Set up for runge Kutta
    yn = np.array((e1,gamma,a1))
    constants = np.array((K,J1,J2,I,mu1,M))
    T_seconds = Tint*365*24*3600
    
    
    #Call RK solver
    output = RungeKutta(yn,constants,T_seconds,fs)
    
    print ('Numerical evolution has completed')
    
    return output


def get_orbital_evolution(m0,m1,f1,e1,beta,m2,e2,I,gamma,Tint,fs):
    
    print ('Getting the orbit')
    
    #Calculate some constants and useful values
    K,J1,J2,mu1,M,a1 = setup(m0,m1,m2,f1,e1,e2,beta)


    #Define time variable
    T_seconds = Tint*365*24*3600
    t = np.arange(0,T_seconds,1/fs)
  
    
    y = np.array((e1,gamma,a1))
    constants = np.array((K,J1,J2,I,mu1,M))
    
    

     #Get approx a

    da = adot(y,constants)
    approx_a = da*t + a1
    

    #Get approx gamma
    dg = gdot(y,constants)
    approx_gamma = t*dg + gamma
    

    
    #Approx e
    ge = e1*(1-e1**2)**(-5/2) * (1 + 121*e1**2/304)
    he = ge/(e1*(1-e1**2))
    psi = get_psi(approx_a,approx_gamma,da,dg,a1,gamma,t)
    AA = 5*K*(1-np.cos(I)**2)/J1
    Cprime = -64*G**3 * mu1 * M**2 / (5 * c**5)

    CC = np.log(e1) - 0.5*np.log(1-e1**2)
    
    alpha_KL = AA*psi + CC
    alpha_GW = -19/36 * Cprime*he/da * (approx_a**(-3) - a1**(-3))

    alpha = alpha_KL + alpha_GW

    
    approx_e = np.exp(alpha)/np.sqrt(1+np.exp(2*alpha))

    
    #output
    out = np.zeros((len(t),4))
    out[:,0] = t
    out[:,1] = approx_e
    out[:,2] = approx_gamma
    out[:,3] = approx_a
    
    
    #out[:,0] = t
    #out[:,1] = e1
    #out[:,2] = gamma
    #out[:,3] = a1
    
    
    return out
# ---------- All below this line are internal functions----------




def derivs(y,constants):

    e_dot = edot(y,constants)
    gamma_dot = gdot(y,constants)
    a_dot = adot(y,constants)


    return np.array((e_dot, gamma_dot,a_dot))



def RungeKutta(yn,constants,Tint,fs):
    
    #Setup timing precision
    h = 1/fs
    t = 0
    nsteps = int(Tint*fs) + 1
    
    
    #Define output array
    out = np.zeros((nsteps,4)) #t,e,gamma,a
    counter = 0
    out[counter,0] = t
    out[counter,1] = yn[0]
    out[counter,2] = yn[1] 
    out[counter,3] = yn[2] 
    counter = counter + 1

    
    print (h, Tint)
    while t < Tint:
        
    

        
        k1 = h * derivs(yn,constants)
        k2 = h * derivs(yn+k1/2,constants)
        k3 = h * derivs(yn+k2/2,constants)
        k4 = h * derivs(yn+k3,constants)

        
        ynew = yn + (k1 + 2*k2 + 2*k3 + k4)/6
        yn = ynew
        
   
    
        t = t + h
        if counter < nsteps:
            out[counter,0] = t
            out[counter,1] = yn[0]
            out[counter,2] = yn[1]  
            out[counter,3] = yn[2] 
        
        
            counter = counter + 1    
    
    return out






def get_psi(a,gamma,da,dg,a0,g0,t):
    
    part1 = -2*a0**2*np.cos(2*gamma)*dg**2
    part2 = 2*a0*da*dg*(-2*t*np.cos(2*gamma)*dg + np.sin(2*gamma))
    part3 = da**2 * (np.cos(2*gamma) - 2*t**2*np.cos(2*gamma)*dg**2 + 2*t*dg*np.sin(2*gamma))
    
    return (part1 + part2+part3)/(4*dg**3)




def setup(m0,m1,m2,f1,e1,e2,beta):
#Setup initial constants, sma etc for hierarchical triple 


    #Masses
    M = m0+m1
    mu1 = (m0*m1)/(m0+m1)

    #Inner Binary calcs
    eps1 = 1-e1**2
    a1 = sma(m0+m1,f1)#semi major axis via K3
    J1 = ang_mom(m0,m1,a1)
    
    #Outer Binary
    eps2 = 1-e2**2
    a2 = beta*a1
    J2 = ang_mom(m0+m1,m2,a2)


    
    #Other
    Kprime = 3*G*m0*m1*m2 / (8*(m0+m1) * a2**3*(1-e2**2)**(3/2))
    

    
    
    return Kprime,J1,J2,mu1,M,a1






#---Keplerian functions for angular momentum and semi major axis
def sma(M,f):
    #Total mass, orbital frequency
    return (G*M/(2*np.pi*f)**2)**(1/3)

def ang_mom(m1,m2,a):
    M = m1+m2
    return m1*m2 * np.sqrt(G*a/M)
