from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from EvalF import Fn
from OrbitalMotion import setup
from OrbitalDerivatives import *
import sys
from scipy.optimize import curve_fit



AU = 1.49e11 #meters




def plot_compare_motion(data1,data2):
    
    fig = plt.figure(figsize=(24,10))
    ax1 = plt.subplot2grid((3,2), (0,0))
    ax2 = plt.subplot2grid((3,2), (1,0), sharex=ax1)
    ax3 = plt.subplot2grid((3,2), (2,0),sharex=ax1)
    
    ax4 = plt.subplot2grid((3,2), (0,1))
    ax5 = plt.subplot2grid((3,2), (1,1),sharex=ax4)
    ax6 = plt.subplot2grid((3,2), (2,1),sharex=ax4)
    
    ax = [ax1,ax2,ax3]
    axR = [ax4,ax5,ax6]
    
    
    
    
    
    data_to_figure(data1,ax,'k')
    data_to_figure(data2,ax,'C0')
    differences(data1,data2,axR)
    

    
    fs = 25
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax4.get_xticklabels(), visible=False)
    plt.setp(ax5.get_xticklabels(), visible=False)


    ax1.tick_params(axis='both', which='major', labelsize=fs)
    ax2.tick_params(axis='both', which='major', labelsize=fs)
    ax3.tick_params(axis='both', which='major', labelsize=fs)
    ax4.tick_params(axis='both', which='major', labelsize=fs)
    ax5.tick_params(axis='both', which='major', labelsize=fs)
    ax6.tick_params(axis='both', which='major', labelsize=fs)


    ax3.set_xlabel('t [years]',fontsize=fs)
    ax6.set_xlabel('t [years]',fontsize=fs)



    ax1.set_ylabel('$e$', fontsize = fs)
    ax2.set_ylabel(r'$\gamma$', fontsize = fs)
    ax3.set_ylabel(r'$a$ [AU]', fontsize = fs)
    

    
    
def data_to_figure(data,ax,c):
    
    t = data[:,0] / (365*24*3600)
    e1 = data[:,1]
    g1 = data[:,2]
    a1 = data[:,3] / AU
    
    
    ax[0].plot(t,e1,c=c)
    
    ax[1].plot(t,g1,c=c)
    
    #if c == 'k':
      #  ax[1].plot(t,g1,c=c)
    #else:
        #ax[1].plot(t,np.sin(2*g1),c=c)
    
    
    
    
    ax[2].plot(t,a1,c=c)

    
    
    
def differences(data,data1,ax):
    
    t = data[:,0]  / (365*24*3600)
    de = (data[:,1] - data1[:,1]) / data1[:,1]
    dg = (data[:,2] - np.sin(2*data1[:,2]))
    da = (data[:,3] - data1[:,3]) / AU
    
    
    ax[0].plot(t,de)
    ax[1].plot(t,dg)
    ax[2].plot(t,da)
    
    
    
    
    
def plot_motion(data1):
    
    
    
    print ('Plotting the orbital parameter evolution')
    
    
    fig = plt.figure(figsize=(24,10))
    ax1 = plt.subplot2grid((3,2), (0,0))
    ax2 = plt.subplot2grid((3,2), (1,0), sharex=ax1)
    ax3 = plt.subplot2grid((3,2), (2,0),sharex=ax1)
    
    ax4 = plt.subplot2grid((3,2), (0,1))
    ax5 = plt.subplot2grid((3,2), (1,1), sharex=ax4)
    ax6 = plt.subplot2grid((3,2), (2,1),sharex=ax4)
    
    

    
    t = data1[:,0] / (365*24*3600)
    e1 = data1[:,1]
    g1 = data1[:,2]
    a1 = data1[:,3] #/ AU
    

    
    maxe = max(e1)
    mine = min(e1)
    ax1.axhline(maxe,linestyle='--')
    ax1.axhline(mine,linestyle='--')
    






    ax1.plot(t,e1)
    ax2.plot(t,g1)
    ax3.plot(t,a1)

    
    ax4.plot(t,e1)
    ax5.plot(t,g1)
    ax6.plot(t,a1)

    
    fs = 25
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax4.get_xticklabels(), visible=False)
    plt.setp(ax5.get_xticklabels(), visible=False)


    ax1.tick_params(axis='both', which='major', labelsize=fs)
    ax2.tick_params(axis='both', which='major', labelsize=fs)
    ax3.tick_params(axis='both', which='major', labelsize=fs)
    ax4.tick_params(axis='both', which='major', labelsize=fs)
    ax5.tick_params(axis='both', which='major', labelsize=fs)
    ax6.tick_params(axis='both', which='major', labelsize=fs)


    ax3.set_xlabel('t [years]',fontsize=fs)
    ax6.set_xlabel('t [years]',fontsize=fs)



    ax1.set_ylabel('$e$', fontsize = fs)
    ax2.set_ylabel(r'$\gamma$', fontsize = fs)
    ax3.set_ylabel(r'$a$ [AU]', fontsize = fs)
        
    
    t_upper = 0.1 #upper limit of t for zooming in
    ax4.set_xlim(-0.01,t_upper)
    
    
    for i in np.arange(len(t)):
        ti = t[i]
        if ti > t_upper:
            ind = i
            break
  

    
    try:
        ax4.set_ylim(min(e1[0:ind]),max(e1[0:ind]))
        ax5.set_ylim(min(g1[0:ind]),max(g1[0:ind]))
        ax6.set_ylim(min(a1[0:ind]),max(a1[0:ind]))
    except:
        pass
        
     
    
    plt.subplots_adjust(hspace=-0.01)
    
    
    
    
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


        
    #f1 = orbital frequency
    t_upper = 5 
    ax3.set_xlim(0,t_upper)
    plt.subplots_adjust(hspace=-0.01)
    

    
    #Save figures
    extent1 = ax1.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    extent2 = ax2.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    
    points = extent1.get_points()
    x0 = points[0][0] * 0.5
    x1 = points[1][0] *1.05
    y1 = points[1][1]
    points = extent2.get_points()
    y0 = points[0][1]
    
    
    my_blit_box = Bbox(np.array([[x0,y0],[x1,y1]]))

    plt.savefig('../../Manuscript/figures/GW_WaveformsLONG.png', dpi=300,bbox_inches=my_blit_box.expanded(1.0, 1.25))

    
    
    ax3.set_ylabel('$h_{+} $', fontsize = fs)
    ax4.set_ylabel(r'$h_{\times}$', fontsize = fs)
    extent1 = ax3.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    extent2 = ax4.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    
    points = extent1.get_points()
    x0 = points[0][0] * 0.9
    x1 = points[1][0] *1.05
    y1 = points[1][1]
    points = extent2.get_points()
    y0 = points[0][1]
    
    
    my_blit_box = Bbox(np.array([[x0,y0],[x1,y1]]))

    plt.savefig('../../Manuscript/figures/GW_WaveformsSHORT.png', dpi=300,bbox_inches=my_blit_box.expanded(1.0, 1.25))
    
    
def plot_GW_frequency(f,h1,h2, S):
    fig = plt.figure(figsize=(14,10))
    ax1 = plt.subplot2grid((1,1), (0,0))
    
    ax1.loglog(f,h1)
    ax1.loglog(f,h2)
    ax1.loglog(f,S)   

    
    fs = 25
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    


    ax1.tick_params(axis='both', which='major', labelsize=fs)

    
    
    
    
    
    
    
    
def func_cos(t,g,om,d,bit):
    return g*om*np.cos(om*t + bit) + d
    
    
def plot_derivatives(data,constants,analytics):
    fig = plt.figure(figsize=(14,10))
    ax1 = plt.subplot2grid((1,1), (0,0))
    
    #Numerical
    t = data[:,0]
    e = data[:,1]
    g = data[:,2]
    a = data[:,3]
    y = np.array((e,g,a))
    dg = gdot(y,constants)
    
    
    omega = 2*dg[0]
    print ('omega =',omega)
    test = np.sin(omega*t + 2*g[0])
    ax1.plot(t,test,'--', c='r')
    
    ga = analytics[:,2]
    
    output=np.zeros((len(g),3)) 
    output[:,0] = t
    output[:,1] = np.sin(2*g)
    output[:,2] = test
    np.savetxt('example.txt',output)
    
    ax1.plot(t,np.sin(2*g),c='C0')
    ax1.plot(t,np.sin(2*ga),c='k')
    
    
    #Analytical
    #e = analytics[:,1]
    ga = analytics[:,2]
    #a = analytics[:,3]
    #y = np.array((e,g,a))
    #dg_approx = gdot(y,constants)
    

    
    #fs = 30
    #ax1.tick_params(axis='both', which='major', labelsize=fs)
    #ax1.plot(t,dg)
    #ax1.plot(t,dg_approx)
    #ax1.plot(t,func_cos(t, *popt))

    
    
    
    
    
    