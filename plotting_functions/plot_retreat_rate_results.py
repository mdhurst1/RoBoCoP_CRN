# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:10:34 2016

Script to plot the results of RockyCoastCRN experiments to explore the influence
of thinning beaches on the concentrations of CRNs built up in the 
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
plt.figure(1,figsize=(6,4))

ERates = [2,5,10,15,20,-99]


for i in range (0,len(ERates)):
    if (ERates[i] != -99):
        FileName = "../driver_files/" + "E" + str(ERates[i]) + "_Ch0_Bh0_Bw0.pdat"
    else:
        FileName = "../driver_files/" + "E5_Ch50_Bh1_Bw50.pdat"
        
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
    
    #setup color
    Color = float(i)/float(len(ERates))
    
    #plot
    if (ERates[i] != -99):
        plt.plot(X,N,'-',color=cm.Paired(Color),lw=2,label=str(ERates[i])+"cm yr$^{-1}$")
    else:
        plt.plot(X,N,'-',color='k',lw=2,label=r"5 cm yr$^{-1}$; $B_h$ = 1 m; $B_w$ = 50 m")

#modify plot
plt.ylabel('$^{10}$Be Concentration (atoms g$^{-1}$)')
plt.xlabel('Distance Offshore (m)')
plt.xlim(-50,1000)
plt.yticks([0,2000,4000,6000,8000])

#Display legend
plt.rcParams.update({'legend.labelspacing':0.1}) 
plt.rcParams.update({'legend.columnspacing':1.0}) 
plt.rcParams.update({'legend.numpoints':1}) 
plt.rcParams.update({'legend.frameon':False}) 
plt.rcParams.update({'legend.handlelength':1.0}) 
plt.rcParams.update({'legend.fontsize':8})
plt.legend(loc=1,ncol=1)

plt.savefig("./RetreatRates.pdf")

plt.show()