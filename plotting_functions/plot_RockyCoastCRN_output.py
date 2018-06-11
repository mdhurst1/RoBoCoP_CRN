# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:10:34 2016

@author: mhurst
"""

#import modules
import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

#open morphology file and read
FileName = "../driver_files/RockyCoastCRN.dat"

# open the file
f = open(FileName,'r')

# get the lines and find out how many lines
Lines = f.readlines()
NoLines = len(Lines)

#Get header info and setup X coord
Header = Lines[0].strip().split(" ")
NXNodes = float(Header[0])
PlatformWidth = float(Header[1])

# setting up the x coordinates
X = np.linspace(0,PlatformWidth+1,NXNodes)
    
#setup figure
fig = plt.figure(1,figsize=(6,6))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

# choose colour map
ColourMap = cm.YlGnBu_r
MaxTime = float(Lines[1].strip().split(" ")[0])

#Loop through data and plot
for i in range(1,NoLines,3):
    
    #Get data    
    ZLine = Lines[i+1].strip().split(" ")
    NLine = Lines[i+2].strip().split(" ")
    Time = float(ZLine[0])
    Z = np.array(ZLine[1:],dtype="float64")
    N = np.array(NLine[1:],dtype="float64")
    
    #mask for NDVs
    mask = Z != -9999
    Zplot = Z[mask]
    Xplot = X[mask]
    Nplot = N[mask]
    
    ax1.plot(Xplot,Zplot,'-',c=ColourMap(Time/MaxTime))
    ax2.plot(Xplot,Nplot,'-',c=ColourMap(Time/MaxTime))
    
ax2.set_xlabel("Distance (m)")
ax2.set_ylabel(r"$^{10}$Be Concentration")
ax1.set_ylabel("Elevation (m)")

plt.savefig("test_output.png", dpi=300)  
