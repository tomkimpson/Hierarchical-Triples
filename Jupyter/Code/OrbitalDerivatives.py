from __future__ import division
import numpy as np
import sys
c = 3e8
G = 6.67e-11
Msolar = 2e30
AU = 1.49e11 #meters

def adot(y,constants):
    
    e = y[0]
    a = y[2]
    
    mu1 = constants[4]
    M = constants[5]
    
    da_GW = -64/5 * G**3*mu1*M**2/(c**5 * a**3) * (1-e**2)**(-7/2) * (1 + 73*e**2 / 24 + 37*e**4 / 96)


    return 0.

    

def edot(y,constants):
    e = y[0]
    gamma = y[1]
    a = y[2]
    
    K = constants[0]
    J1 = constants[1]
    I = constants[3]
    mu1 = constants[4]
    M = constants[5]
    
    de_Kozai = 5*K*e*(1-e**2)*a**2*(1-np.cos(I)**2)*np.sin(2*gamma)/J1 
    de_GW = -304/15 * G**3*mu1*M**2/(c**5*a**4) * e/(1-e**2)**(5/2) * (1+121*e**2/304)
    
    #print ('Big A a^2 =', 5*K*a**2*(1-np.cos(I)**2)/J1, e)
    #sys.exit()


    return de_Kozai #+ de_GW


def edot_linear(y,constants):
    #Only the oscillatory Kozai bit of edot
    e = y[0]
    gamma = y[1]
    a = y[2]
    
    K = constants[0]
    J1 = constants[1]
    I = constants[3]
    mu1 = constants[4]
    M = constants[5]
    
    de_Kozai = 5*K*e*(1-e**2)*a**2*(1-np.cos(I)**2)*np.sin(2*gamma)/J1 
    de_GW = -304/15 * G**3*mu1*M**2/(c**5*a**4) * e/(1-e**2)**(5/2) * (1+121*e**2/304)

    
    
    return 0


def gdot(y,constants):
    
    e = y[0]
    gamma = y[1]
    a = y[2]
    
    K = constants[0]
    J1 = constants[1]
    J2 = constants[2]
    I = constants[3]
    mu1 = constants[4]
    M = constants[5]
    
    
    dg_part1 = 2*(1-e**2) - 5*(1-e**2 - np.cos(I)**2) *np.sin(gamma)**2
    dg_part2 = (1 - e**2 + 5*e**2 * np.cos(gamma)**2)*np.cos(I)
    dg_PN=3/(c**2*a*(1-e**2)) * (G*M/a)**(3/2)
    
    
    dg_KL = 2*K*a**2*(dg_part1/J1 + dg_part2/J2)

    
    '''
    u = 1-e**2
    DA1 = 4*K*a**2*u/J1
    DA2 = -10*K*a**2/J1 *u*np.sin(gamma)**2
    DA3 = 10*K*a**2/J1 *np.cos(I)**2*np.sin(gamma)**2
    
    DB1 = 2*K*a**2*np.cos(I)/J2 * u
    DB2 = 10*K*a**2*np.cos(I)/J2 * np.sin(gamma)**2
    DB3 = -10*K*a**2*np.cos(I)/J2 * u*np.sin(gamma)**2
    
    decomp = DA1 + DA2 + DA3 + DB1+DB2+DB3
    '''
    
    
    
    
    return dg_KL + dg_PN
