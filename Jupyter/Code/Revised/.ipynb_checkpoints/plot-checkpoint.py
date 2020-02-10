from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from scipy.optimize import curve_fit   


c = 3e8
G = 6.67e-11
Msolar = 2e30
AU = 1.49e11 #meters
pc = 3.086e16 #parsec in m


def plot_motion_and_derivatives(motion,const):
    
    
    
    print ('Plotting the orbital parameter evolution')
    
    
    fig = plt.figure(figsize=(24,20))
    ax1 = plt.subplot2grid((4,2), (0,0))
    ax2 = plt.subplot2grid((4,2), (1,0), sharex=ax1)
    ax3 = plt.subplot2grid((4,2), (2,0),sharex=ax1)
    ax4 = plt.subplot2grid((4,2), (3,0), sharex=ax1)
    
    ax5 = plt.subplot2grid((4,2), (0,1))
    ax6 = plt.subplot2grid((4,2), (1,1), sharex=ax5)
    ax7 = plt.subplot2grid((4,2), (2,1),sharex=ax5)
    #ax8 = plt.subplot2grid((4,2), (3,1), sharex=ax5)
    
    

    #Motion
    t = motion[:,0] / (365*24*3600)
    e1 = motion[:,1]
    g1 = motion[:,2]
    a1 = motion[:,3] #/ data1[0,3]
    J1 = motion[:,4] #/ data1[0,4]
    
    #Derivatives
    #-edot
    u = 1-e1**2
    C = const[2]
    A = const[3]
    eKL = A*e1*u*a1**2 * np.sin(2*g1)/J1
    eGW = 19/12 * C * a1**(-4) * e1 * u**(-5/2)*(1+121/304 * e1**2)
    edot = eKL + eGW
    
    #-gdot
    K = const[0]
    J2 = const[1]
    C = const[2]
    A = const[3]
    I = const[5]
    Lambda = const[6]
    
    part1 = 2*u - 5*(u-np.cos(I)**2)*np.sin(g1)**2
    part2 = np.cos(I)*(u + 5*e1**2 * np.sin(g1)**2)
    
    gKL = 2*K*a1**2*(part1/J1 + part2/J2)
    gPN = Lambda * a1**(-2.5)/u
    
    
    gdot = gKL + gPN
    
    #-adot
    C = const[2]
    fe = 1 + 73/24 * e1**2 + 37/96 * e1**4  
    adot = C*a1**(-3) *u**(-7/2) * fe
    
    #-Jdot
    eta = const[4]
    he = 1 + 7/8 * e1**2
    Jdot = eta*a1**(-7/2)*u**(-2) * he
 
    
    
    
    #also get the orbital frequency - assumes total mass is 60 solar masses
    f1 = np.sqrt(G*60*Msolar *a1**(-3) / (4*np.pi**2))
    
    
    
    #Plot it all
    
    
    
    maxe = max(e1)
    mine = min(e1)
    ax1.axhline(maxe,linestyle='--')
    ax1.axhline(mine,linestyle='--')
    print ('E limits =', maxe,mine, maxe-mine)

    ax1.plot(t,e1)
    ax2.plot(t,np.sin(2*g1))
    ax3.plot(t,a1/a1[0])
    ax4.plot(t,f1)
    #ax4.plot(t,J1/J1[0])
    
    ax5.plot(t,edot)
    ax6.plot(t,gdot)
    ax7.plot(t,adot)
    #ax8.plot(t,Jdot)


    
    fs = 25
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax5.get_xticklabels(), visible=False)
    plt.setp(ax6.get_xticklabels(), visible=False)
    plt.setp(ax7.get_xticklabels(), visible=False)


    ax1.tick_params(axis='both', which='major', labelsize=fs)
    ax2.tick_params(axis='both', which='major', labelsize=fs)
    ax3.tick_params(axis='both', which='major', labelsize=fs)
    ax4.tick_params(axis='both', which='major', labelsize=fs)
    
    
    ax5.tick_params(axis='both', which='major', labelsize=fs)
    ax6.tick_params(axis='both', which='major', labelsize=fs)
    ax7.tick_params(axis='both', which='major', labelsize=fs)
   # ax8.tick_params(axis='both', which='major', labelsize=fs)
   #
   

    ax4.set_xlabel('t [years]',fontsize=fs)




    ax1.set_ylabel('$e$', fontsize = fs)
    ax2.set_ylabel(r'$\sin( 2 \gamma)$', fontsize = fs)
    ax3.set_ylabel(r'$a$ [AU]', fontsize = fs)
    ax4.set_ylabel(r'$f$ [Hz]', fontsize = fs)
        

    ax5.set_ylim(-0.00001,0.000015)
    ax3.set_ylim(0,1.1)
    
    ax4.set_yscale('log')
    
    #ax1.set_xlim(-0.001, 0.1)
    #ax5.set_xlim(-0.001, 0.1)

    ax6.set_ylim(0.00001,0.00005)
    
    plt.subplots_adjust(hspace=-0.01)











def plot_motion(data1):
    
    
    
    print ('Plotting the orbital parameter evolution')
    
    
    fig = plt.figure(figsize=(14,10))
    ax1 = plt.subplot2grid((4,1), (0,0))
    ax2 = plt.subplot2grid((4,1), (1,0), sharex=ax1)
    ax3 = plt.subplot2grid((4,1), (2,0),sharex=ax1)
    ax4 = plt.subplot2grid((4,1), (3,0), sharex=ax1)
    
    
    

    
    t = data1[:,0] / (365*24*3600)
    e1 = data1[:,1]
    g1 = data1[:,2]
    a1 = data1[:,3] / data1[0,3]
    J1 = data1[:,4] / data1[0,4]
    
    
    maxe = max(e1)
    mine = min(e1)
    ax1.axhline(maxe,linestyle='--')
    ax1.axhline(mine,linestyle='--')
    print ('E limits =', maxe,mine, maxe-mine)

    ax1.plot(t,e1)
    ax2.plot(t,np.sin(2*g1))
    ax3.plot(t,a1)
    ax4.plot(t,J1)


    
    fs = 25
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax3.get_xticklabels(), visible=False)

    ax1.tick_params(axis='both', which='major', labelsize=fs)
    ax2.tick_params(axis='both', which='major', labelsize=fs)
    ax3.tick_params(axis='both', which='major', labelsize=fs)
    ax4.tick_params(axis='both', which='major', labelsize=fs)
   

    ax4.set_xlabel('t [years]',fontsize=fs)




    ax1.set_ylabel('$e$', fontsize = fs)
    ax2.set_ylabel(r'$\sin( 2 \gamma)$', fontsize = fs)
    ax3.set_ylabel(r'$a$ [AU]', fontsize = fs)
    ax4.set_ylabel(r'$J_1$', fontsize = fs)
        

    #ax1.set_xlim(-0.01,0.2)
    
    plt.subplots_adjust(hspace=-0.01)
    
    
def plot_compare_motion(data1,data2):
    
    fig = plt.figure(figsize=(24,10))
    ax1 = plt.subplot2grid((3,2), (0,0))
    ax2 = plt.subplot2grid((3,2), (1,0), sharex=ax1)
    ax3 = plt.subplot2grid((3,2), (2,0), sharex=ax1)

    
    ax5 = plt.subplot2grid((3,2), (0,1))
    ax6 = plt.subplot2grid((3,2), (1,1))
    ax7 = plt.subplot2grid((3,2), (2,1))

    
    
    ax = [ax1,ax2,ax3]
    axR = [ax5,ax6,ax7]
    
    
    
    
    
    data_to_figure(data1,ax,'k')
    data_to_figure(data2,ax,'C0')
    differences(data1,data2,axR)
    

    
    fs = 25
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
   
    
    plt.setp(ax5.get_xticklabels(), visible=False)
    plt.setp(ax6.get_xticklabels(), visible=False)
 


    ax1.tick_params(axis='both', which='major', labelsize=fs)
    ax2.tick_params(axis='both', which='major', labelsize=fs)
    ax3.tick_params(axis='both', which='major', labelsize=fs)
    
    ax5.tick_params(axis='both', which='major', labelsize=fs)
    ax6.tick_params(axis='both', which='major', labelsize=fs)
    ax7.tick_params(axis='both', which='major', labelsize=fs)


    ax3.ticklabel_format(useOffset=False)
    ax3.set_xlabel('t [years]',fontsize=fs)
    ax7.set_xlabel('t [years]',fontsize=fs)



    ax1.set_ylabel('$e$', fontsize = fs)
    ax2.set_ylabel(r'$\sin (2\gamma)$', fontsize = fs)
    ax3.set_ylabel(r'$a$ [mAU]', fontsize = fs)
    
    print ('LIM CHANGES')
    lim = 0.025
    ax5.set_ylim(-lim,lim)
    ax6.set_ylim(-lim,lim)
    ax7.set_ylim(-lim,lim)
    
    path = '/Users/tomkimpson/PhD/PI/PI_Work/Manuscript/figures/'
    plt.savefig(path+'compare_canonical.png',dpi=300)
        
    

    
    
    
def data_to_figure(data,ax,c):
    
    t = data[:,0] / (365*24*3600)
    e1 = data[:,1]
    g1 = data[:,2]
    a1 = data[:,3] / AU
    a1 =a1 * 1e3 #milli AU

    
    
    
    ax[0].plot(t,e1,c=c)
    ax[1].plot(t,np.sin(2*g1),c=c)  
    ax[2].plot(t,a1,c=c)

        
    
    
    
    
    
def differences(data,data1,ax):
    
    t = data[:,0]  / (365*24*3600)
    
    de = (data[:,1] - data1[:,1]) / data1[:,1]
    
    #dg = (np.sin(data[:,2]) - np.sin(data1[:,2])) #/ np.sin(data1[:,2])
    dg = (data[:,2] - data1[:,2]) / data1[:,2]
    
    print (data[:,2])
    print (data1[:,2])
    print (dg)
    
          
    da = (data[:,3] - data1[:,3]) / data1[:,3]
    
    print ('The average percentage error in a was:', np.mean(da)*100)
    
          
          
 
    
    
    ax[0].plot(t,de)
    ax[1].plot(t,dg)
    ax[2].plot(t,da)

    
    
    
    
def temp(data):
    
    #Plot environment setup
    fig = plt.figure(figsize=(24,20))
    ax1 = plt.subplot2grid((2,1), (0,0))
    ax2 = plt.subplot2grid((2,1), (1,0),sharex=ax1)
    
    #Load the data
    t = data[:,0] 
    e = data[:,1]
    g1 = data[:,2]
    a = data[:,3] 
    J1 = data[:,4] 
    
    
    #Plot the derivatives of gamma
    u = 1 - e**2
    K =  1.1565942657526629e+24
    lam =  2.3869282539699428e+16
    I = 1.0471975511965976
    J2 =  7.182511894390675e+46


    P1 = 2*u - 5*np.sin(g1)**2*(u-np.cos(I)**2)
    P2 = (u + 5*e**2*np.sin(g1)**2)*np.cos(I)
    box = P1/J1 + P2/J2
    #g0 = 2*K*a**2*box+ lam * a**(-2.5) * u**(-1) #the form before we rearrange it
    
    NSP = 2*K*a**2 * u * (2/J1 + np.cos(I)/J2) + lam * a**(-2.5) * u**(-1)
    SP = 10*K*a**2*(e**2*np.cos(I)/J2 - (u-np.cos(I)**2)/J1)
    
    gALT = NSP + SP*np.sin(g1)**2
    
    ax1.plot(t,gALT, linestyle='--', label ='analytical')
    
    
    
    #And get the numerical derivative just out of interest
    dx = t[-1]-t[-2]
    dt = []
    for i in range(len(t)-1):
        inter = (t[i+1] - t[i])/2 + t[i]
        dt.extend([inter])
    
 
    dy = np.diff(g1)/dx
    ax1.plot(dt,dy, label='Numerical diff')
    
    
    #Plot gamma
    ax2.plot(t, g1, c='g', label='true')
      
        
    #Now, your simple model uses a linear model for gamma where the omega,c = 2.17209725e-05 6.85614205e-01
    omega, c = 2.17209725e-05, 6.85614205e-01
    gamma_approx = omega*t + c
    
    ax1.plot(t,omega*np.ones(len(t)), label='linear')
    ax2.plot(t,gamma_approx, label='linear solution')

    
    
    #Now can we fit the derivative?
    approx_derivative = extract(t,gALT)
    ax1.plot(t,approx_derivative,label = 'approximation')
    
    #And how does the integral of the approximate derivative do?

    
    B_A, B_omega, B_offset, B_D,gamma_f, gamma_g,H = 3.14193996e-05,  4.33566266e-05, -1.86768534e-01, -3.70351335e-05, 2.16994938e-05, 6.86571691e-01, 3.23403676e-05
    
    integration_constant = g1[0] - integral(B_A, B_omega, B_offset, B_D,gamma_f, gamma_g,H,0) 
    approx_gamma_oscil = integral(B_A, B_omega, B_offset, B_D,gamma_f, gamma_g,H,t) + integration_constant
    ax2.plot(t,approx_gamma_oscil,label = 'oscil solution')
    
    
    
    ## - end
    ax1.legend()
    ax2.legend()
    
    
    
#-------------------



def integral(A,omega,C,D,f,g,H,t):
    
    part1 = -A*np.cos(C+t*(omega-2*f)-2*g)/(4*(2*f-omega))
    
    part2 =  A*np.cos(C+t*(omega+2*f)+2*g)/(4*(2*f+omega))
    
    part3 = 0.5*A*(np.sin(C)*np.sin(t*omega)/omega - np.cos(C)*np.cos(omega*t)/omega)
    
    part4 = D*(f*t + g)/(2*f)
    
    part5 = -D*np.sin(2*(f*t + g))/(4*f)
    
    part6 = H*t
    
   
    return part1 + part2 +part3 + part4 + part5 + part6



def doube_trig_function(t,B_A, B_omega, B_offset, B_D,gamma_f, gamma_g,H):
    B = B_A*np.sin(B_omega*t + B_offset) +B_D
    s2g = np.sin(gamma_f*t + gamma_g)**2
    return H + B*s2g
    
def extract(t,f):
    print ('Extract for double trig func')
    func = doube_trig_function
    
    
    H= 3.631335880308346e-05
    B_A = 2.50885879e-05
    B_omega = 4.32315349e-05
    B_offset = -1.76843846e-01 
    B_D = -2.17565784e-05
    gamma_f = 2.17209725e-05
    gamma_g = 6.85614205e-01
    
    
    
    
    
    p0 = (B_A, B_omega, B_offset, B_D,gamma_f, gamma_g,H)
        
    popt, pcov = curve_fit(func, t,f,p0=p0)
    print ('G approx:', popt)
    return func(t, *popt)





#---------------------
    
    
    
    
    
    
    
    
    
def compare_GW(data1,data2):    
    
    
    fig = plt.figure(figsize=(24,10))
    ax1 = plt.subplot2grid((2,2), (0,0))
    ax2 = plt.subplot2grid((2,2), (1,0), sharex=ax1)
    
   
    
    
    tyear = data1[:,0] / (365*24*3600)
    hp1 = data1[:,1]
    hc1 = data1[:,2]
    hp2 = data2[:,1]
    hc2 = data2[:,2]
    
    dh_plus = abs((hp1 - hp2)/hp1)
    dh_cross = abs((hc1 - hc2)/hc1)
    
    ax1.plot(tyear,dh_plus, c='C0')
    ax2.plot(tyear,dh_plus,c='C1')

    

    
    fs = 25
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    




    ax1.tick_params(axis='both', which='major', labelsize=fs)
    ax2.tick_params(axis='both', which='major', labelsize=fs)
    


    ax2.set_xlabel('t [years]',fontsize=fs)

    ax1.set_ylabel('$dh_{+} $', fontsize = fs)
    ax2.set_ylabel('$dh_{\times} $', fontsize = fs)
    
    ax1.set_yscale('log')
    ax2.set_yscale('log')

    ax1.set_xlim(0.09,0.11)
    
    
def plot_GW(data1,f1):
    
    print ('Plotting the GW')
    
    fig = plt.figure(figsize=(24,10))
    ax1 = plt.subplot2grid((2,2), (0,0))
    ax2 = plt.subplot2grid((2,2), (1,0), sharex=ax1)
    
    ax3 = plt.subplot2grid((2,2), (0,1))
    ax4 = plt.subplot2grid((2,2), (1,1),sharex=ax3)
    
    
  
    
    tyear = data1[:,0] / (365*24*3600)
    torb = data1[:,0] * f1
    hplus = data1[:,1]
    hcross = data1[:,2]
    
    
    print ('Max h plus =', max(hplus) )
    ax1.axhline(1.2623854965230435e-21, linestyle='--')
    
    
    ax1.plot(tyear,hplus, c='C0')
    ax2.plot(tyear,hcross,c='C1')

    ax3.plot(torb,hplus,c='C0')
    ax4.plot(torb,hcross,c='C1')
    

    
    fs = 25
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax3.get_xticklabels(), visible=False)



    ax1.tick_params(axis='both', which='major', labelsize=fs)
    ax2.tick_params(axis='both', which='major', labelsize=fs)
    ax3.tick_params(axis='both', which='major', labelsize=fs)
    ax4.tick_params(axis='both', which='major', labelsize=fs)
   


    ax2.set_xlabel('t [years]',fontsize=fs)
    ax4.set_xlabel(r'$N_{orbits}$',fontsize=fs)



    ax1.set_ylabel('$h_{+} $', fontsize = fs)
    ax2.set_ylabel(r'$h_{\times}$', fontsize = fs)


    
    
    ax1.set_xlim(0.09,0.11)
        
    #f1 = orbital frequency
    t_upper = 5 
    ax3.set_xlim(0,t_upper)
    plt.subplots_adjust(hspace=-0.01)
    
    
    #path = '/Users/tomkimpson/PhD/PI/PI_Work/Manuscript/figures/'
    #plt.savefig(path+'GW_canonical.png',dpi=300)
    

    
    #Save figures
    #extent1 = ax1.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    #extent2 = ax2.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    
    #points = extent1.get_points()
    #x0 = points[0][0] * 0.5
    #x1 = points[1][0] *1.05
    #y1 = points[1][1]
    #points = extent2.get_points()
    #y0 = points[0][1]
    
    
    #my_blit_box = Bbox(np.array([[x0,y0],[x1,y1]]))

    #plt.savefig('../../Manuscript/figures/GW_WaveformsLONG.png', dpi=300,bbox_inches=my_blit_box.expanded(1.0, 1.25))

    
    
    #ax3.set_ylabel('$h_{+} $', fontsize = fs)
    #ax4.set_ylabel(r'$h_{\times}$', fontsize = fs)
    #extent1 = ax3.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    #extent2 = ax4.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    
    #points = extent1.get_points()
    #x0 = points[0][0] * 0.9
    #x1 = points[1][0] *1.05
    #y1 = points[1][1]
    #points = extent2.get_points()
    #y0 = points[0][1]
    
    
    #my_blit_box = Bbox(np.array([[x0,y0],[x1,y1]]))

    #plt.savefig('../../Manuscript/figures/GW_WaveformsSHORT.png', dpi=300,bbox_inches=my_blit_box.expanded(1.0, 1.25))
    
    
    
def plot_GW_frequency(f,h1,h2, S):
    fig = plt.figure(figsize=(15,15))
    ax1 = plt.subplot2grid((1,1), (0,0))
    
    ax1.loglog(f,h1, 'C2')
    ax1.loglog(f,h2,'C1')
    ax1.loglog(f,S, C='C0')   

    
    fs = 25
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    


    ax1.tick_params(axis='both', which='major', labelsize=fs)
    ax1.set_ylabel('$h(f) $', fontsize = fs)
    ax1.set_xlabel(r'f [Hz]', fontsize = fs)
    ax1.tick_params(axis='both', which='major', labelsize=fs)
    
    path = '/Users/tomkimpson/PhD/PI/PI_Work/Manuscript/figures/'
    plt.savefig(path+'GW_overlap.png',dpi=300)

    
    
    
    

    
def plot_compare_motionSPLIT(numerical,sec1,sec2):
    
    fig = plt.figure(figsize=(24,10))
    
    ax1 = plt.subplot2grid((1,2), (0,0))
    ax5 = plt.subplot2grid((1,2), (0,1))
    year = (365*24*3600)
    
    #Plot the full numerical solution
    ax1.plot(numerical[:,0]/year,numerical[:,1])
    
    #Plot section1
    ax1.plot(sec1[:,0]/year,sec1[:,1])
    
    #Plot section2
    ax1.plot(sec2[:,0]/year,sec2[:,1])
    
    
    #Get the relative differences
    
    #section 1
    lim = int(len(numerical)/2)
    xN = numerical[0:lim+1,0]
    yN = numerical[0:lim+1,1]
    de = (yN - sec1[:,1]) / yN
    ax5.plot(xN, de)
    
    #section 2
    xN = numerical[lim:-1,0]
    yN = numerical[lim:-1,1]
    de = (yN - sec2[:,1]) / yN
    ax5.plot(xN, de)

    
        

    
    
  
    #formatting
    fs = 25
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    

    ax1.tick_params(axis='both', which='major', labelsize=fs)
    ax5.tick_params(axis='both', which='major', labelsize=fs)
   


