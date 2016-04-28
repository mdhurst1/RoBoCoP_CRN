# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:10:34 2016

Script to plot the results of RockyCoastCRN experiments to explore the influence
of beaches on the concentrations of CRNs built up in the 
platform surface

Martin Hurst,
Feb 10th 2016

@author: mhurst
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rc

# Customise figure style #
rc('font',size=8)
rc('ytick.major',pad=5)
rc('xtick.major',pad=5)
padding = 5
plt.figure(1,figsize=(6,8))
ax1 = plt.axes([0.1,0.775,0.85,0.2])
ax2 = plt.axes([0.1,0.48,0.85,0.2])
ax3 = plt.axes([0.1,0.1,0.85,0.35])

FileNames = ["BeachWidth_Bh0_Bw0","BeachWidth_Bh1_Bw50","Variable_Bw50_Bh1"]
Labels = ["$B_H$ = 0 m; $B_W$ = 0 m;","$B_H$ = 1 m; $B_W$ = 50 m;","$B_H$ = 1 m; $B_W$ = 50 +/- 30 m;"]

Time = np.arange(0.,1000.,1.)
BeachWidth = 50.+30.*np.sin(2.*np.pi*Time/100.)
ax1.plot(Time,BeachWidth,'k-')

for i in range (0,len(FileNames)):
    FileName = "../driver_files/" + FileNames[i] + ".pdat"
    f = open(FileName,'r')
    Lines = f.readlines()
    NoLines = len(Lines)
    
    #Get header info and setup X coord
    Header = Lines[0].strip().split(" ")
    NXNodes = float(Header[0])
    PlatformWidth = float(Header[1])
    X = np.linspace(0,PlatformWidth+1,NXNodes)
    
    Z = np.array(Lines[-3].strip().split(" "), dtype="float64")[1:]
    N = np.array(Lines[-1].strip().split(" "), dtype="float64")[1:]
    
    #plt.plot(X,Z,'k-')
    
    #plot timeseries of platform evolution
    for j in range(10,NoLines-1,8):
            
        #Get data    
        Line = Lines[j].strip().split(" ")
        Time = float(Line[0])
        ZRock = np.array(Line[1:],dtype="float64")
        
        #mask for NDVs
        mask = ZRock != -9999
        Zplot = ZRock[mask]
        Xplot = X[mask]
        
        if (i == 0):
            ax2.plot([Xplot[0]-20,Xplot[0],Xplot[0]],[5,5,Zplot[0]],'k-')
            ax2.plot(Xplot,Zplot,'k-',zorder=10)
            ax2.text(Xplot[0],Zplot[0]+4.5, str(Time/1000.) + " ka",rotation=-90)
    
        #plot variable beach widths
        else:
            #define platform surface
            ZPlatform = -(Xplot-Xplot[0])*0.01
            #get profile of minimum beach width
            ZBeachMin = Zplot
            ZBeachMin[(Xplot-Xplot[0]) < 20] = 1.
            ZBeachMin[(Xplot-Xplot[0]) > 20] = 1.-0.125*np.power((Xplot[(Xplot-Xplot[0]) > 20]-Xplot[0]-20.),2./3.)
            mask = ZBeachMin > ZPlatform
            XBeachMin = Xplot[mask]
            ZBeachMin = ZBeachMin[mask]
            #get profile of maximum beach width
            ZBeachMax = Zplot
            ZBeachMax[(Xplot-Xplot[0]) < 80] = 1.
            ZBeachMax[(Xplot-Xplot[0]) > 80] = 1.-0.125*np.power((Xplot[(Xplot-Xplot[0]) > 80]-Xplot[0]-80.),2./3.)
            mask = ZBeachMax > ZPlatform
            XBeachMax = Xplot[mask]
            ZBeachMax = ZBeachMax[mask]
            
            #ax2.plot(Xplot,Zplot,'k-')
            ax2.plot(XBeachMin,ZBeachMin,'k-')
            ax2.plot(XBeachMax,ZBeachMax,'k-')
            
            XFill = np.concatenate((XBeachMax,XBeachMax[::-1]))
            ZFill = np.concatenate((ZBeachMax,ZPlatform[mask][::-1]))
            ax2.fill(XFill,ZFill,color=[1.0,1.0,0.5])
            
    if (i == 0):
        ax3.plot(X,N,'k-',lw=2,label=str(Labels[i]),zorder=10)
    else:
        ax3.plot(X,N,'k-',lw=2,label=str(Labels[i]))

ax1.set_ylabel('Beach Width $B_w$ (m)')
ax1.set_xlabel('Time (Years)')
ax1.set_ylim(0,100)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')
ax1.text(40,90,'(a)')

ax2.set_ylabel('Elevation (m)')
ax2.set_xlim(-50,1000)
ax2.set_ylim(-10,8)
ax2.set_xticklabels([])
ax2.set_yticks([-10,-5,0,5])
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.yaxis.set_ticks_position('left')
ax2.set_xticklabels([])
ax2.set_xticks([])
ax2.text(0,6.5,'(b)')


ax3.set_xlim(-50,1000)
ax3.set_yticks([0,250,500,750,1000,1250,1500,1750])
ax3.set_ylabel('$^{10}$Be Concentration (atoms g$^{-1}$)')
ax3.set_xlabel('Distance Offshore (m)')
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.yaxis.set_ticks_position('left')
ax3.xaxis.set_ticks_position('bottom')
ax3.text(0,1700,'(c)')

##Display legend
#plt.rcParams.update({'legend.labelspacing':0.1}) 
#plt.rcParams.update({'legend.columnspacing':1.0}) 
#plt.rcParams.update({'legend.numpoints':1}) 
#plt.rcParams.update({'legend.frameon':False}) 
#plt.rcParams.update({'legend.handlelength':0.5}) 
#plt.legend(loc=1,ncol=1,title="Beach Parameters")
#leg = plt.gca().get_legend()
#
##set fontsize to small
#ltext  = leg.get_texts()
#plt.setp(ltext, fontsize=8) 
#
plt.xlim(-50,1000)

plt.savefig("./VariableBeachWidth.pdf")

plt.show()