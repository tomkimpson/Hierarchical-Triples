from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
AU = 1.49e11 #meters



def plot_motion(data1):
    
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
    a1 = data1[:,3] / AU

    
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
  

    ax4.set_ylim(min(e1[0:ind]),max(e1[0:ind]))
    ax5.set_ylim(min(g1[0:ind]),max(g1[0:ind]))
    ax6.set_ylim(min(a1[0:ind]),max(a1[0:ind]))
    
    
    plt.subplots_adjust(hspace=-0.01)
    
    
    
    
def plot_GW(data1,f1):
    fig = plt.figure(figsize=(24,10))
    ax1 = plt.subplot2grid((2,2), (0,0))
    ax2 = plt.subplot2grid((2,2), (1,0), sharex=ax1)
    
    ax3 = plt.subplot2grid((2,2), (0,1))
    ax4 = plt.subplot2grid((2,2), (1,1),sharex=ax3)
    
    
    
    
    t = data1[:,0] / (365*24*3600)
    hplus = data1[:,1]
    hcross = data1[:,2]
    
    ax1.plot(t,hplus, c='C0')
    ax2.plot(t,hcross,c='C1')

    ax3.plot(t,hplus,c='C0')
    ax4.plot(t,hcross,c='C1')

    
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
    ax4.set_xlabel('t [years]',fontsize=fs)



    ax1.set_ylabel('$h_{+} $', fontsize = fs)
    ax2.set_ylabel(r'$h_{\times}$', fontsize = fs)


        
    #f1 = orbital frequency
    Torbit = 1/f1 / (365*24*3600)
    t_upper = 5 * Torbit #upper limit of t for zooming in
    ax3.set_xlim(0,t_upper)
    
    
    
    
    plt.subplots_adjust(hspace=-0.01)
    
    
    
    
def plot_GW_frequency(f,hf, S, hphen, fphen):
    fig = plt.figure(figsize=(14,10))
    ax1 = plt.subplot2grid((1,1), (0,0))
    
    ax1.loglog(f,hf)
    ax1.loglog(f,S)   
    ax1.loglog(fphen,hphen)
    fs = 25
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    


    ax1.tick_params(axis='both', which='major', labelsize=fs)
