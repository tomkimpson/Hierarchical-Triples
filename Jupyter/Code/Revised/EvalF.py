from __future__ import division 
from numpy import sin as Sin 
from numpy import pi as Pi 
from numpy import sqrt as Sqrt 
from numpy import arctan as ArcTan 
from numpy import cos as Cos 
from numpy import tan as Tan 
from numpy import angle as Arg 

#This expression was evaluated via a Taylor exampsion up to order n= 4  
                                                                       
#about the point e0 = 0.5  
                          
def Fn(A,B,omega,DD,t,e0): 
    output1 =  (64.8884812316277*omega*t + 942.0732843034996*A**2*omega*t + 629.5038475185829*A**4*omega*t + 942.0732843034996*B**2*omega*t + 1259.0076950371658*A**2*B**2*omega*t + 629.5038475185829*B**4*omega*t - 561.9549219909673*DD*omega*t - 4259.767899697493*A**2*DD*omega*t - 4259.767899697493*B**2*DD*omega*t + 1884.1465686069992*DD**2*omega*t + 5036.030780148663*A**2*DD**2*omega*t + 5036.030780148663*B**2*DD**2*omega*t - 2839.845266464995*DD**3*omega*t + 1678.676926716221*DD**4*omega*t + A*(561.9549219909673 + A**2*(2129.8839498487464 - 5036.030780148663*DD) + B**2*(2129.8839498487464 - 5036.030780148663*DD) - 3768.2931372139983*DD + 8519.535799394986*DD**2 - 6714.707706864884*DD**3)*Cos(omega*t) + A*B*(-942.0732843034996 - 839.3384633581105*A**2 - 839.3384633581105*B**2 + 4259.767899697493*DD - 5036.030780148663*DD**2)*Cos(2.*omega*t) - 236.65377220541623*A**3*Cos(3.*omega*t) + 709.961316616249*A*B**2*Cos(3.*omega*t) + 559.5589755720737*A**3*DD*Cos(3.*omega*t) - 1678.676926716221*A*B**2*DD*Cos(3.*omega*t) + 209.83461583952763*A**3*B*Cos(4.*omega*t) - 209.83461583952763*A*B**3*Cos(4.*omega*t) - 561.9549219909673*B*Sin(omega*t) - 2129.883949848747*A**2*B*Sin(omega*t) - 2129.883949848747*B**3*Sin(omega*t) + 3768.293137213998*B*DD*Sin(omega*t) + 5036.030780148664*A**2*B*DD*Sin(omega*t) + 5036.030780148664*B**3*DD*Sin(omega*t) - 8519.535799394987*B*DD**2*Sin(omega*t) + 6714.707706864884*B*DD**3*Sin(omega*t) - 471.03664215174973*A**2*Sin(2.*omega*t) - 419.66923167905526*A**4*Sin(2.*omega*t) + 471.03664215174973*B**2*Sin(2.*omega*t) + 419.66923167905526*B**4*Sin(2.*omega*t) + 2129.883949848747*A**2*DD*Sin(2.*omega*t) - 2129.883949848747*B**2*DD*Sin(2.*omega*t) - 2518.015390074332*A**2*DD**2*Sin(2.*omega*t) + 2518.015390074332*B**2*DD**2*Sin(2.*omega*t) + 709.961316616249*A**2*B*Sin(3.*omega*t) - 236.65377220541632*B**3*Sin(3.*omega*t) - 1678.676926716221*A**2*B*DD*Sin(3.*omega*t) + 559.5589755720739*B**3*DD*Sin(3.*omega*t) + 52.45865395988191*A**4*Sin(4.*omega*t) - 314.7519237592915*A**2*B**2*Sin(4.*omega*t) + 52.45865395988191*B**4*Sin(4.*omega*t))/omega
    return output1