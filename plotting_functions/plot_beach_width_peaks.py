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
plt.figure(1,figsize=(3.3,3.0))

FileNames = ["BeachWidth_Bh0_Bw0","BeachWidth_Bh1_Bw0","BeachWidth_Bh1_Bw10","BeachWidth_Bh1_Bw20","BeachWidth_Bh1_Bw50","BeachWidth_Bh1_Bw100"]
Labels = ["$B_H$ = 0 m; $B_W$ = 0 m;", "$B_H$ = 1 m; $B_W$ = 0 m;","$B_H$ = 1 m; $B_W$ = 10 m;","$B_H$ = 1 m; $B_W$ = 20 m;","$B_H$ = 1 m; $B_W$ = 50 m;","$B_H$ = 1 m; $B_W$ = 100 m;"]
BeachWidths = np.array([0.,10.,20.,50.,100.])
Ns = np.zeros(len(BeachWidths))

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
    
    Ind = np.argmax(N)
    if (i>0):
        Ns[i-1] = N[Ind]

#linear regression
A = np.vstack([BeachWidths, np.ones(len(BeachWidths))]).T
m, c = np.linalg.lstsq(A,Ns)[0]

ModelN = m*BeachWidths+c

plt.plot(BeachWidths,Ns,'ko')
plt.plot(BeachWidths,ModelN,'k--')

plt.xlabel('Beach Width $B_w$ (m)')
plt.ylabel('Maximum $^{10}$Be Concentration (atoms g$^{-1}$)')
plt.savefig("../BeachWidthPeaks.pdf")
plt.xlim(-5,105)
plt.ylim(1000,1800)
plt.tight_layout()
plt.show()