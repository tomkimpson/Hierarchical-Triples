from __future__ import division
import numpy as np
from scipy import special as sp
import sys
G = 6.67e-11
c=3e8
def GW(orbital_motion,constants):
    
    print ('Getting the waveform')
    t = orbital_motion[:,0]
    e = orbital_motion[:,1]
    g = orbital_motion[:,2]
    a = orbital_motion[:,3]
    
    M = constants[0] #total mass of inner binary
    nmodes = constants[1]
    iota = constants[2]
    mu = constants[3]
    D = constants[4]
    
    
    MA =  np.sqrt(G*M) * a**(-3/2)*t #mean anomaly
    fdynamic = np.sqrt(G*M/(4*np.pi**2)) * a**(-3/2) #reparam of a to f
    omega = 2*np.pi*fdynamic
    
    #print (fdynamic)
    #print ('First:', (M * G / c**2 * np.pi * fdynamic[0])**(2/3))
    #sys.exit()
    
    #P = mu/D  #prefactor of mu / D
    #P = 1
    #AA = P * (2*np.pi*fdynamic*M)**(2/3) 
    
    #AA = mu/D
    
    
    #Convert to geometric units
    mu = mu * G/c**2
    M = M * G/c**2
    D = D
    omega = omega/c
    
    AA = mu/D * (M*omega)**(2/3)
    #AA = mu/D
    #AA = (M*omega)**(2/3)
        
        
    #print ('AA = ', AA)  
    
    
    
    
    
    
    
    
    
    
    waveform_out = np.zeros((len(t), 3))
    waveform_out[:,0] = t 
    

    
    
    
    for n in np.arange(1,nmodes+1):
        J_2 = sp.jv(n-2,n*e)
        J_1 = sp.jv(n-1,n*e) 
        Jn = sp.jv(n,n*e) 
        J1 = sp.jv(n+1,n*e)
        J2 = sp.jv(n+2,n*e)
        
        an = -n*AA*(J_2 - 2*e*J_1 + 2*Jn/n + 2*e*J1 - J2)*np.cos(n*MA)
        bn = -n*AA*np.sqrt((1-e**2)) * (J_2 - 2*Jn + J2)*np.sin(n*MA)
        cn = 2*AA*Jn*np.cos(n*MA)
                
            
        hplus = -(1+np.cos(iota)) * (an*np.cos(2*g) - bn*np.sin(2*g)) + (1-np.cos(iota)**2)*cn
        hcross = 2*np.cos(iota)*(bn*np.cos(2*g) + an*np.sin(2*g))
        
        waveform_out[:,1] = waveform_out[:,1] + hplus
        waveform_out[:,2] = waveform_out[:,2] + hcross
        

        
     
    return waveform_out
        
        