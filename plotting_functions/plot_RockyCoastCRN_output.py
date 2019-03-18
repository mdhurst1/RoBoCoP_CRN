# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:10:34 2016

@author: mhurst
"""

#import modules
import matplotlib
#matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

#setup figure
fig = plt.figure(1,figsize=(6,6))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)


# choose colour map
ColourMap = cm.hot

# get filename
FileName = "../driver_files/scalby_8cm_testX.dat"

# open the file
f = open(FileName,'r')

# get the lines and find out how many lines
Lines = f.readlines()
NoLines = len(Lines)

#Get header info and setup X coord
Header = Lines[0].strip().split(" ")
NXNodes = float(Header[0])
PlatformWidth = float(Header[1])

MaxTime = float(Lines[1].strip().split(" ")[0])

for i in range(1,NoLines,4):
    
    #Get data  
    X = np.array(Lines[i].strip().split(" ")[1:], dtype="float64")
    ZPlatform = np.array(Lines[i+1].strip().split(" ")[1:], dtype="float64")
    ZBeach = np.array(Lines[i+2].strip().split(" ")[1:], dtype="float64")
    N = np.array(Lines[i+3].strip().split(" ")[1:], dtype="float64")
    Time = float(Lines[i].strip().split(" ")[0])
    
    #mask for NDVs
    mask = ZBeach != -9999
    Zbeach = ZBeach[mask]
    Zplat = ZPlatform[mask]
    Xplot = X[mask]
    Nplot = N[mask]

    Colour = Time/MaxTime
    ax1.plot(Xplot,Zbeach,'--',c=ColourMap(Colour))
    ax1.plot(Xplot,Zplat,'-',c=ColourMap(Colour))
    ax2.plot(Xplot,Nplot,'-',c=ColourMap(Colour))


ax2.set_xlabel("Distance (m)")
ax2.set_ylabel(r"$^{10}$Be Concentration")
ax1.set_ylabel("Elevation (m)")

#limit x axis to 250m
#plt.xlim(0,250) 

plt.savefig("test_output_SY_8cm_testX.png", dpi=300)  
plt.show()