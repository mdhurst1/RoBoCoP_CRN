# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:10:34 2016

Script to plot the results of RoBoCoP and RockyCoastCRN

Martin Hurst,
March 7th 2016

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

#First plot the morphology through time
FileName = "../driver_files/ShoreProfile.xz"
f = open(FileName,'r')
Lines = f.readlines()
NoLines = len(Lines)

ax1 = plt.subplot(211)

#Get header info and setup X coord
Header = Lines[0].strip().split(" ")
Times = []
RetreatRate = []
for j in range(1,NoLines,2):
    
    Line = (Lines[j].strip().split(" "))
    Time = float(Line[0])
    X = np.array(Line[1:],dtype="float64")
    Z = np.array((Lines[j+1].strip().split(" "))[1:],dtype="float64")
    
    ax1.plot(X,Z,'k-',lw=6)

ax1.set_xticklabels([])
plt.ylabel("Elevation (m)")

#now plot CRN concentration through time
FileName = "../driver_files/CRNConcentrations.xn"
f = open(FileName,'r')
Lines = f.readlines()
NoLines = len(Lines)

ax2 = plt.subplot(212)
#Get header info and setup X coord
Header = Lines[0].strip().split(" ")
for j in range(1,NoLines,2):
    X = np.array((Lines[j].strip().split(" "))[1:],dtype="float64")
    N = np.array((Lines[j+1].strip().split(" "))[1:],dtype="float64")
    
    ax2.plot(X,N,'r-',lw=2)
    
plt.xlabel('Distance (m)')
plt.ylabel('Concentration (atoms g$^{-1}$)')

plt.tight_layout()
plt.show()