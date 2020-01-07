from __future__ import division
import numpy as np
from scipy import special as sp
import scipy.integrate




G = 6.67e-11
c=3e8

def Gwaves(motion,constants):
    
    print ('Getting the waveform')
    t = motion[:,0]
    e = motion[:,1]
    g = motion[:,2]
    a = motion[:,3]
    
    M = constants[0] #total mass of inner binary
    nmodes = constants[1]
    iota = constants[2]
    mu = constants[3]
    D = constants[4]
    
    
    MA =  np.sqrt(G*M) * a**(-3/2)*t #mean anomaly
    fdynamic = np.sqrt(G*M/(4*np.pi**2)) * a**(-3/2) #reparam of a to f
    omega = 2*np.pi*fdynamic
    

    
    #Convert to geometric units
    mu = mu * G/c**2
    M = M * G/c**2
    D = D
    omega = omega/c
    
    AA = mu/D * (M*omega)**(2/3)

        
    
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




def overlap(data1,data2):
    
    #Get data in Fourier regime
    f1,hp1,hc1 = FT(data1)
    f2,hp2,hc2 = FT(data2)
  
    f = f1
    
    #Noise
    S = noise(f)
    
    
    
    hplusN = hp1
    hcrossN = hc1
    
    hN = np.sqrt(abs(hplusN)**2 + abs(hcrossN)**2) #numerical
    hsig1 = hN
    normN = 2 * scipy.integrate.simps( (hN * np.conj(hN) + np.conj(hN) * hN)/(S),f)
    hN = hN / np.sqrt(normN)
    #overlap = 2 * scipy.integrate.simps( (hN * np.conj(hN) + np.conj(hN) * hN)/(S),f)


    hplusA = hp2
    hcrossA = hc2
    
    hA = np.sqrt(abs(hplusA)**2 + abs(hcrossA)**2) #numerical
    hsig2 = hA
    normA = 2 * scipy.integrate.simps( (hA * np.conj(hA) + np.conj(hA) * hA)/(S),f)
    hA = hA / np.sqrt(normA)
    
    
    overlap = 2 * scipy.integrate.simps( (hN * np.conj(hA) + np.conj(hN) * hA)/(S),f)

 


    
    print ('overlap = ', overlap)
    
    
    return f,hsig1,hsig2, np.sqrt(S)

    
    
    
    
def FT(data):    
    
    t = data[:,0]
    hplus = data[:,1]
    hcross = data[:,2]
    
    
    dt = t[1] - t[0]
    fs = 1/dt
    print (fs)

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
    
    return f, hp, hc
    
    
    
    
    
    
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