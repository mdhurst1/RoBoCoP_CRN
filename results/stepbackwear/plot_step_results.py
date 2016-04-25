# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:10:34 2016

Script to plot the results of RockyCoastCRN experiments to explore the influence
of platform erosion processes on the concentrations of CRNs built up in the 
platform surface

Martin Hurst,
Feb 9th 2016

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

#setup figure
plt.figure(1,figsize=(6,8),facecolor="white")

#Filename lists
FileNames = ["Steps100cm","Steps50cm","Steps20cm","Steps10cm","NoSteps"]
StepSizes = [100,50,20,10,0]

#setup subplot for plotting  concentrations
ax1 = plt.axes([0.1,0.1,0.8,0.32])
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')
plt.xlabel("Distance (m)")
plt.ylabel("$^{10}$Be Concentration (atoms g$^{-1}$ yr$^{-1}$)")
for i in range (0,len(FileNames)):
    
    #setup subplot
    ax = plt.axes([0.1,0.88-float(i)*0.11,0.8,0.09])
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.set_yticks([-10,-5,0,5])
    
    #setup color
    Color = float(i)/float(len(FileNames))
    
    FileName = "../data/blockremoval/BlockRemoval_" + FileNames[i] + ".pdat"
    f = open(FileName,'r')
    Lines = f.readlines()
    NoLines = len(Lines)
    
    #Get header info and setup X coord
    Header = Lines[0].strip().split(" ")
    NXNodes = float(Header[0])
    PlatformWidth = float(Header[1])
    X = np.linspace(0,PlatformWidth+1,NXNodes)
    
    for j in range(10,NoLines,8):
        #Get data    
        Line = Lines[j].strip().split(" ")
        Time = float(Line[0])
        Z = np.array(Line[1:],dtype="float64")
        
        #mask for NDVs
        mask = Z != -9999
        Zplot = Z[mask]
        Xplot = X[mask]
        
        plt.plot([Xplot[0]-20,Xplot[0],Xplot[0]],[5,5,0],'-',color=cm.Paired(Color))
        if (i == 0):
            plt.text(Xplot[0]-20,Zplot[0]+8, str(Time/1000.) + " ka")
        plt.plot(Xplot,Zplot,'-',color=cm.Paired(Color))

    #Offset -= 15
    if (i ==0):
        ax.text(950,5,'(a)')
    if (i == 2):
        plt.ylabel('Elevation (m)') 
    plt.xlim(-50,1000)

    #Plot concentrations
    #Get data    
    Line = Lines[NoLines-1].strip().split(" ")
    Time = float(Line[0])
    N = np.array(Line[1:],dtype="float64")
    #plot concentration data    
    ax1.plot(X,N,'-',color=cm.Paired(Color),linewidth=2,label=str(StepSizes[i]/100.) + " m")


    
#Display legend
plt.rcParams.update({'legend.labelspacing':0.1}) 
plt.rcParams.update({'legend.columnspacing':1.0}) 
plt.rcParams.update({'legend.numpoints':1}) 
plt.rcParams.update({'legend.frameon':False}) 
plt.rcParams.update({'legend.handlelength':0.5})
plt.rcParams.update({'legend.fontsize':8})
ax1.legend(loc=1,ncol=1,title="Bed Spacing")
ax1.text(0,2800,'(b)')
ax1.set_xlim(-50,1000)
ax.set_xticklabels([])

plt.savefig("../figures/BlockRemoval.pdf")
plt.show()