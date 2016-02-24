# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 10:22:56 2016

Calculate Time spent at each elevation in a tidal cycle

Martin Hurst, 24/2/16

@author: mhurst
"""
import numpy as np
import matplotlib.pyplot as plt

#Tides
H = 2.5
dh = 0.1
T=12.
Time = np.arange(0.,T,0.001)
WaterLevels = H*np.sin((2.*np.pi*Time/(T)));
NTidalValues = len(WaterLevels)

#bin the tidal elevations
ElevationBins = np.arange(-H,H+dh,dh)
hist = np.histogram(WaterLevels,bins=ElevationBins,density=True)
WaterLevels2 = hist[1]
Weights = hist[0]

#convert to erosion distribution
#Calculate integrated erosion for a tidal cycle
Z = np.arange(-10,10,0.1)
Ez = np.zeros(len(Z))
for i in range(0,len(WaterLevels2)-1):
#for i in range(0,20):
    #Calculate the surf force (kg/m^2) reaching the water line
    SurfForce = Weights[i]
        
    #Calculate erosion across the platform
    ErosionMask = Z <= WaterLevels2[i]
    WaterDepth = WaterLevels2[i]-Z[ErosionMask]
    
    #Do erosion here including exponential decline with water depth
    Ez[ErosionMask] += SurfForce*np.exp(-WaterDepth)

plt.figure(1)
plt.plot(Time,WaterLevels)

#plt.plot(t,Z,'k-',lw=2)

plt.figure(2)
plt.barh(ElevationBins[:-1], hist[0]/10.,height=dh)

plt.figure(3)
plt.plot(Ez,Z,'k-')
plt.show()
