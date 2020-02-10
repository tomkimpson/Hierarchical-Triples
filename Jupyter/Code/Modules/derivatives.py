from __future__ import division
import numpy as np

#----- Turn on/off variou effects
KL = 1
PN = 1
GW = 1



#-------
c = 3e8
G = 6.67e-11
Msolar = 2e30
AU = 1.49e11 #meters
pc = 3.086e16 #parsec in m
    

def edot(y,const):
    e = y[0]
    g = y[1]
    a = y[2]
    J1 = y[3]
    
    u = 1-e**2
    
    C = const[2]
    A = const[3]

    
    eKL = A*e*u*a**2 * np.sin(2*g)/J1
    eGW = 19/12 * C * a**(-4) * e * u**(-5/2)*(1+121/304 * e**2)
    
    
    if KL == 0:
        eKL = 0.0
    if GW == 0:
        eGW = 0.0
    
    return eKL + eGW
    
    
def gdot(y,const):
    
    e = y[0]
    g = y[1]
    a = y[2]
    J1 = y[3]
    
    u = 1-e**2
    
    
    K = const[0]
    J2 = const[1]
    C = const[2]
    A = const[3]
    I = const[5]
    Lambda = const[6]
    
    part1 = 2*u - 5*(u-np.cos(I)**2)*np.sin(g)**2
    part2 = np.cos(I)*(u + 5*e**2 * np.sin(g)**2)
    
    
    
    gKL = 2*K*a**2*(part1/J1 + part2/J2)
    gPN = Lambda * a**(-2.5)/u
    
    
   
 
    
    

    
    if KL == 0:
        gKL = 0.0
    if PN == 0:
        gPN = 0.0
        
    
    return gKL + gPN




def adot(y,const):
    
    e = y[0]
    a = y[2]
    u = 1-e**2

    C = const[2]
    
    fe = 1 + 73/24 * e**2 + 37/96 * e**4  
    
    aGW = C*a**(-3) *u**(-7/2) * fe
    
    if GW == 0:
        aGW = 0.0
    
    return aGW



def Jdot(y,const):
    
    e = y[0]
    a = y[2]
    u = 1-e**2
    
    eta = const[4]
    he = 1 + 7/8 * e**2
    
    JGW = eta*a**(-7/2)*u**(-2) * he
    
    if GW == 0:
        JGW = 0.0
    
    return JGW