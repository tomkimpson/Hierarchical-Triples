from __future__ import division 
from numpy import sin as Sin 
from numpy import cos as Cos 

#This expression was evaluated a Taylor exampsion up to order n= 2  
                                                                   
def Fn(ebar,edot,e0,omega,t): 
    
    print (ebar)
    print (omega)
    print (Cos(t))
    #bit1 = 342.625*t + 1109.9583333333333*edot*t**2 + 1664.7777777777776*t*(3*ebar**2 + 2*edot**2*t**2)
    bit2= - (2219.9166666666665*ebar*Cos(omega*t))/omega 
    #bit3 = -(19977.333333333332*ebar*edot*t*Cos(omega*t))/omega 
    #bit4 = (19977.333333333332*ebar*edot*Sin(omega*t))/omega**2  
    #bit5 = -(2497.1666666666665*ebar**2*Sin(2*omega*t))/omega
    
    #output1 = bit1+bit2+bit3+bit4+bit5
    output1=bit2
    return output1