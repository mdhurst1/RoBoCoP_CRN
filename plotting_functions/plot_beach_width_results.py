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
plt.figure(1,figsize=(6,6))
ax1 = plt.axes([0.1,0.65,0.85,0.3])
ax2 = plt.axes([0.1,0.1,0.85,0.5])
FileNames = ["BeachWidth_Bh0_Bw0","BeachWidth_Bh1_Bw0","BeachWidth_Bh1_Bw10","BeachWidth_Bh1_Bw20","BeachWidth_Bh1_Bw50","BeachWidth_Bh1_Bw100"]
Labels = ["$B_H$ = 0 m; $B_W$ = 0 m;", "$B_H$ = 1 m; $B_W$ = 0 m;","$B_H$ = 1 m; $B_W$ = 10 m;","$B_H$ = 1 m; $B_W$ = 20 m;","$B_H$ = 1 m; $B_W$ = 50 m;","$B_H$ = 1 m; $B_W$ = 100 m;"]

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
        ZRockplot = ZRock[mask]
        Xplot = X[mask]
        
        #setup color
        Color = float(i)/float(len(FileNames))
        ax1.plot(Xplot,ZRockplot,'-',color=cm.Paired(Color))
        
        
        if (i == 0):
            ax1.plot([Xplot[0]-20,Xplot[0],Xplot[0]],[5,5,ZRockplot[0]],'k-')
            ax1.plot(Xplot,ZRockplot,'k-',zorder=10)
            ax1.text(Xplot[0],ZRockplot[0]+4, str(Time/1000.) + " ka",rotation=-90)
        
        
    if (i == 0):
        ax2.plot(X,N,'k-',lw=2,label=str(Labels[i]),zorder=10)
    else:
        ax2.plot(X,N,'-',color=cm.Paired(Color),lw=2,label=str(Labels[i]))

ax1.set_ylabel('Elevation (m)')
ax1.set_xlim(-50,1000)
ax1.set_ylim(-10,8)
ax1.set_xticklabels([])
ax1.set_yticks([-10,-5,0,5])
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.yaxis.set_ticks_position('left')
ax1.set_xticklabels([])
ax1.set_xticks([])
ax1.text(0,6.5,'(a)')

ax2.set_xlim(-50,1000)
ax2.set_yticks([0,250,500,750,1000,1250,1500,1750])
ax2.set_ylabel('$^{10}$Be Concentration (atoms g$^{-1}$)')
ax2.set_xlabel('Distance Offshore (m)')
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.yaxis.set_ticks_position('left')
ax2.xaxis.set_ticks_position('bottom')
ax2.text(0,1700,'(b)')

#Display legend
plt.rcParams.update({'legend.labelspacing':0.1}) 
plt.rcParams.update({'legend.columnspacing':1.0}) 
plt.rcParams.update({'legend.numpoints':1}) 
plt.rcParams.update({'legend.frameon':False}) 
plt.rcParams.update({'legend.handlelength':0.5}) 
plt.legend(loc=1,ncol=1,title="Beach Parameters")
leg = plt.gca().get_legend()

#set fontsize to small
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=8) 


plt.xlim(-50,1000)

plt.savefig("../BeachWidth.pdf")

plt.show()