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
ColourMap = cm.copper

# get filename
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

MaxTime = float(Lines[1].strip().split(" ")[0])

for i in range(1,NoLines,3):
    
    #Get data  
    X = np.array(Lines[i].strip().split(" ")[1:], dtype="float64")
    Z = np.array(Lines[i+1].strip().split(" ")[1:], dtype="float64")
    N = np.array(Lines[i+2].strip().split(" ")[1:], dtype="float64")
    Time = float(Lines[i].strip().split(" ")[0])
    
    #mask for NDVs
    mask = Z != -9999
    Zplot = Z[mask]
    Xplot = X[mask]
    Nplot = N[mask]

    Colour = Time/MaxTime
    ax1.plot(Xplot,Zplot,'-',c=ColourMap(Colour))
    ax2.plot(Xplot,Nplot,'-',c=ColourMap(Colour))

ax2.set_xlabel("Distance (m)")
ax2.set_ylabel(r"$^{10}$Be Concentration")
ax1.set_ylabel("Elevation (m)")

plt.savefig("test_output_1.png", dpi=300)  
plt.show()