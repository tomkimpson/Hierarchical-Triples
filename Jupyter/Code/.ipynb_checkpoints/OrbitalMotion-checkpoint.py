from __future__ import division
import numpy as np
from OrbitalDerivatives import *
from EvalF import Fn, edotINTEGRAL
from numpy import sin as Sin
from numpy import cos as Cos


c = 3e8
G = 6.67e-11
Msolar = 2e30
AU = 1.49e11 #meters



def get_orbital_evolutionV2(m0,m1,f1,e1,beta,m2,e2,I,gamma,Tint,fs,numerical_data):
    
    print ('Getting the orbit')

    
    
    #Calculate some constants and useful values
    K,J1,J2,mu1,M,a1 = setup(m0,m1,m2,f1,e1,e2,beta)
    

    


    #Define time variable
    T_seconds = Tint*365*24*3600
    t = np.arange(0,T_seconds,1/fs)
        
    
    y = np.array((e1,gamma,a1))
    constants = np.array((K,J1,J2,I,mu1,M))
    
    #Get some info from the numerical data
    AmpE,omega,C,D = extract(numerical_data,'osc',1)
    
    A1, om1, C1, A2, om2, C2, D = extract3(numerical_data,1)
    eALT = A1*np.sin(om1*t + C1) +A2*np.cos(om2*t + C2) + D

    
    

    

    #A1, om1, C1, A2, om2, C2,A3, om3, C3, D = extract3(numerical_data,2)
    #gALT = A1*np.sin(om1*t + C1) +A2*np.cos(om2*t + C2)+ A3*np.cos(om3*t + C3)+ D

    
    A1, om1, C1, A2, om2, C2,A3, om3, C3, D = extract3(numerical_data,2)
    gALT = A1*np.sin(om1*t + C1) +A2*np.cos(om2*t + C2)+ D

    
    
    #Can we fit gamma?
    OmG, DPHI = extract2(numerical_data[:,0],np.sin(2*numerical_data[:,2]),simple_trig)
    approx_gamma=(OmG*t + DPHI)/2
    
    
    #-------------17/12/2019-----------START
    alpha = 10*K*a1**2*np.cos(I)*(1/J1 + 1/J2)
    CC = -1/np.tan(gamma)
    
    
    LinearGamma = np.arctan(1/(-alpha*t - CC))
    

    
    
    
    
     #-------------17/12/2019-----------END
    
    
   

    
    
    
    
    
    
    
    #Get some derivative info based on initial conditions
    dg = gdot(y,constants)
    de = edot_linear(y,constants)
    
    #Approximate e(t)
    approx_e = AmpE * np.sin(omega*t + C) + D + de*t
    
    #Approximate da
    MidE = max(numerical_data[:,1])- AmpE #Midline of the sinusoid - need this since the oscillations are not symmetric about the initial value
    Cprime = -64*G**3 * mu1 * M**2 / (5 * c**5)
    FN0 = Fn(AmpE,de,MidE,omega,0,D,C)
    FNT = Fn(AmpE,de,MidE,omega,t,D,C)
    Da = a1**4 / 4 - Cprime*FN0
    approx_a = (4*(Cprime*FNT + Da))**(1/4)
    
    #We could also fit a line to the jagged motion as:
    adot, a0 = extract(numerical_data,'lin',3)
    approx_a = adot*t + a0
    approx_a = a1 #No GW emmision

    
    BigA = 5*K*(1-np.cos(I)**2)/J1
    DD = e1 - edotINTEGRAL(AmpE,omega,OmG,DPHI,adot,a0,BigA,0)
    
    newE = edotINTEGRAL(AmpE,omega,OmG,DPHI,adot,a0,BigA,t) + DD
    
    
    
    

    #Get approx gamma, again using a parametric fit
    gd, g0 = extract(numerical_data,'lin',2)
    #gradient = gd*(1+0.1*np.cos(eALT))
    
    y = np.array((0.8,gamma,a1))
    dgBIGE = gdot(y,constants)
    y = np.array((0.35,gamma,a1))
    dgSMALLE = gdot(y,constants)
    
    #print ('dg = ', dgBIGE,dgSMALLE)
    #sys.exit()
    approx_gamma = dg*t + gamma

    


    #output
    out = np.zeros((len(t),4))
    out[:,0] = t
    out[:,1] = eALT
    out[:,2] = approx_gamma
    out[:,3] = approx_a
    

    
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
    
    trange = np.arange(0,Tint,h)
    tfinal = trange[-1]    
    
    t = 0
    nsteps = int(tfinal*fs) + 1

    
    
    
    #Define output array
    out = np.zeros((nsteps,4)) #t,e,gamma,a
    counter = 0
    out[counter,0] = t
    out[counter,1] = yn[0]
    out[counter,2] = yn[1] 
    out[counter,3] = yn[2] 
    counter = counter + 1

    
    
    while t < tfinal:
        
    

        
        k1 = h * derivs(yn,constants)
        k2 = h * derivs(yn+k1/2,constants)
        k3 = h * derivs(yn+k2/2,constants)
        k4 = h * derivs(yn+k3,constants)

        
        ynew = yn + (k1 + 2*k2 + 2*k3 + k4)/6
        yn = ynew
        
   
    
        t = t + h
        
        if counter <= nsteps-1:
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




def doube_trig_function(t,A1,omega1,C1,A2,omega2,C2,D):
    return A1*np.sin(omega1*t + C1) +A2*np.cos(omega2*t + C2) + D

def triple_trig_function(t,A1,omega1,C1,A2,omega2,C2,A3, omega3,C3,D):
    return A1*np.sin(omega1*t + C1) +A2*np.cos(omega2*t + C2)+A3*np.cos(omega3*t + C3) + D


def trig_function(t,A,omega,C,D):
    return A*np.sin(omega*t + C) + D

def linear_function(t,A,B):
    return A*t + B

def simple_trig(t,omega,dphi):
    return np.sin(omega*t+dphi)

def simple_trigHIGHER(t,omega,dphi,omega2):
    return np.sin(omega*t+dphi+omega2*t**2)


def simple_trig2(t,omega,dphi,omega2, dphi2):
    return np.sin(omega*t+dphi) + np.cos(omega2 + dphi2)


def LinOsc(t,omega,dphi,Aprime,omegaprime):
    return np.sin(omega*t+dphi + Aprime*np.sin(omegaprime*t))




from scipy.optimize import curve_fit

def extract(data, method,index):
    t = data[:,0]
    f = data[:,index] 
    
    if method == 'osc':
        
        func = trig_function
        offset = (f.max() + f.min()) / 2
        y_shifted = f - offset
        p0 = (
        (f.max() - f.min()) / 2,
        np.pi * np.sum(y_shifted[:-1] * y_shifted[1:] < 0) / (t.max() - t.min()),
        0,
        offset
        )
    if method == 'lin':
        func = linear_function
        p0 = (1e-6,f[0])

    



    #do a scipy fit
    popt, pcov = curve_fit(func, t,f,p0=p0)
    return popt


def extract2(x,f,func):
    
    offset = (f.max() + f.min()) / 2
    y_shifted = f - offset
    p0 = (
        np.pi * np.sum(y_shifted[:-1] * y_shifted[1:] < 0) / (x.max() - x.min()),
        0
        )
        
    popt, pcov = curve_fit(func, x,f,p0=p0)

    return popt

def extract22(x,f):
    func = LinOsc
    
    offset = (f.max() + f.min()) / 2
    y_shifted = f - offset
    p0 = (
        np.pi * np.sum(y_shifted[:-1] * y_shifted[1:] < 0) / (x.max() - x.min()),
        0,
        0.1,
        1e-4
        )
        
    popt, pcov = curve_fit(func, x,f,p0=p0)

    return popt






def extractHIGHER(x,f,func):
    
    offset = (f.max() + f.min()) / 2
    y_shifted = f - offset
    p0 = (
        np.pi * np.sum(y_shifted[:-1] * y_shifted[1:] < 0) / (x.max() - x.min()),
        0,
        np.sqrt(np.pi * np.sum(y_shifted[:-1] * y_shifted[1:] < 0) / (x.max() - x.min()))
        )
        
    popt, pcov = curve_fit(func, x,f,p0=p0)

    return popt


def extract3(data,index):
    t = data[:,0]
    f = data[:,index] 
    
    
    
    if index == 1:
        func = doube_trig_function
        offset = (f.max() + f.min()) / 2
        y_shifted = f - offset
        p0 = (
        (f.max() - f.min()) / 2,
        np.pi * np.sum(y_shifted[:-1] * y_shifted[1:] < 0) / (t.max() - t.min()),
        0,
        (f.max() - f.min()) / 2,
        np.pi * np.sum(y_shifted[:-1] * y_shifted[1:] < 0) / (t.max() - t.min()),
        0,
        offset
        )
    
    
    if index == 2:
        f = np.sin(2*f)
        func = triple_trig_function
        offset = (f.max() + f.min()) / 2
        y_shifted = f - offset
        p0 = (
        (f.max() - f.min()) / 2,
        np.pi * np.sum(y_shifted[:-1] * y_shifted[1:] < 0) / (t.max() - t.min()),
        0,
        (f.max() - f.min()) / 2,
        np.pi * np.sum(y_shifted[:-1] * y_shifted[1:] < 0) / (t.max() - t.min()),
        0,
        (f.max() - f.min()) / 2,
        np.pi * np.sum(y_shifted[:-1] * y_shifted[1:] < 0) / (t.max() - t.min()),
        0,
        offset
        )

    
    popt, pcov = curve_fit(func, t,f,p0=p0)
    return popt
    
    
    
    
        
def extract4(x,f):
    
    func = LinOsc
    offset = (f.max() + f.min()) / 2
    y_shifted = f - offset
    p0 = (1e-6,f[0],
        (max(np.sin(2*f)) - min(np.sin(2*f))) / 2,
        1e-6,
        0
        )
        
    popt, pcov = curve_fit(func, x,f,p0=p0)
    return popt
    