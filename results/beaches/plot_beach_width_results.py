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
ax = plt.subplot(211)

FileNames = ["Beach_Bh1_Bw0","Beach_Bh1_Bw10","Beach_Bh1_Bw20","Beach_Bh1_Bw50","Beach_Bh1_Bw100","Beach_Bh1_Bw200"]
Labels = ["$B_H$ = 1 m; $B_W$ = 0 m;","$B_H$ = 1 m; $B_W$ = 10 m;","$B_H$ = 1 m; $B_W$ = 20 m;","$B_H$ = 1 m; $B_W$ = 50 m;","$B_H$ = 1 m; $B_W$ = 100 m;","$B_H$ = 1 m; $B_W$ = 200 m;"]
for i in range (0,len(FileNames)):
    FileName = FileNames[i] + "_Morphology.txt"
    f = open(FileName,'r')
    Lines = f.readlines()
    NoLines = len(Lines)
    
    #Get header info and setup X coord
    Header = Lines[0].strip().split(" ")
    NXNodes = float(Header[0])
    PlatformWidth = float(Header[1])
    X = np.linspace(0,PlatformWidth+1,NXNodes)
    
    if i == 0:
        for j in range(5,NoLines-2,4):
            #Get data    
            Line = Lines[j].strip().split(" ")
            Time = float(Line[0])
            ZRock = np.array(Line[1:],dtype="float64")
            
            #mask for NDVs
            mask = ZRock != -9999
            ZRockplot = ZRock[mask]
            Xplot = X[mask]
            
            plt.plot([Xplot[0]-20,Xplot[0],Xplot[0]],[10,10,ZRockplot[0]],'k-')
            if (i == 0):
                plt.text(Xplot[0],ZRockplot[0]+5, str(Time/1000.) + " ka",rotation=-90)
            plt.plot(Xplot,ZRockplot,'k-')
        
        
    Line = Lines[-2].strip().split(" ")
    Time = float(Line[0])
    ZRock = np.array(Line[1:],dtype="float64")
    
    Line = Lines[-1].strip().split(" ")
    Time = float(Line[0])
    ZBeach = np.array(Line[1:],dtype="float64")
    
    #mask for NDVs
    mask = ZRock != -9999
    ZRockplot = ZRock[mask]
    Xplot = X[mask]
    ZBeachplot = ZBeach[mask]
    Color = float(i)/float(len(FileNames))
    plt.plot(Xplot,ZBeachplot,'-',color=cm.Paired(Color))
    
plt.plot([Xplot[0]-20,Xplot[0],Xplot[0]],[10,10,ZRockplot[0]],'k-')
plt.text(Xplot[0],ZRockplot[0]+5, str(Time/1000.) + " ka",rotation=-90)
plt.plot(Xplot,ZRockplot,'k-')
            
plt.ylabel('Elevation (m)')
plt.xlim(-50,1000)
ax.set_xticklabels([])
    
plt.subplot(212)
Color=0
for i in range(0,len(FileNames)):  
    #open morphology file and read
    FileName = FileNames[i] + "_CRNs.txt"
    X, N = np.loadtxt(FileName,skiprows=1, unpack=True)
    Color = float(i)/float(len(FileNames))
    plt.plot(X,N,'-',color=cm.Paired(Color),linewidth=2,label=str(Labels[i]))

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

plt.xlabel('Distance (m)')
plt.ylabel('$^{10}$Be (atoms g$^{-1}$)')
plt.xlim(-50,1000)
plt.tight_layout()
plt.savefig("Beach_Width.pdf")
plt.show()