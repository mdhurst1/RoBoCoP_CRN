# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 14:21:24 2016

Model for the evolution of a shore platform developed for use with
Cosmogenic Radionuclide predictions. Model follows Trenhaile (2000).

Martin D. Hurst
British Geological Survey
February 11th 2016

@author: mhurst
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

#Define constants
M = 6.5*10.**(-3.)          #Conversion coefficient (m^3/kg)
C = 0.5                     #Coeffcient representing breaker type
k = 0.02                    #Coefficient of wave energy attenuation
rho_w = 1025.0              #density of water (kg/m^3)
g = 9.81                    #gravity (m/s^2)
u = 10.0                    #wind velocity (m/s)
MinSurfForce = 100.         #critical resisting force to overcome (kg/m^2)
W = 240.                    #waves/hour!?

#Geometric Parameters
PlatformGradient = 1./10.
#setup intial conditions & spatial parameters
InitialCliffPosition = 500.
InitialCliffHeight = 10.
InitialJunctionElevation = 1.
dz = 0.1
Z = np.arange(-10.,10.,dz)
NoNodes = len(Z)
X = -Z/PlatformGradient+InitialCliffPosition+InitialJunctionElevation/PlatformGradient
X[Z>InitialJunctionElevation] = InitialCliffPosition
X[-1] = 0.

#Tides
TidalAmplitude = 1.0
TidalPeriod=12.42
TidalTimes = np.arange(0,TidalPeriod/2.,0.001)
WaterLevels = -TidalAmplitude*np.cos((2.*np.pi*TidalTimes/(TidalPeriod)));
NTidalValues = len(WaterLevels)

#Sea level
SLR = 0.0001 # m/yr
SeaLevel = 0

#Waves
#significant wave height and period calculated following Kraal et al 2006
#equations 11 and 12, after Pond and Pickard, 1983
Hs = 0.2412*u**2/g      #Significant wave height (m)
T = 7.6924*(u/g)       #peak period (s)

#Incoming wave height (Hb) is related to far field wave height (hs)
#Komar and Gaughan (1972)
Hb = 0.39*(g**0.2)*(T*Hs**2)**0.4

#This can be related to water depth at the point of wave breaking
hb = Hb*1.0/0.78        #breaking wave depth

#setup time control
MaxTime = 10000.
Time = MaxTime
EndTime = 0.
dt = 1.
PlotTime = MaxTime
PlotInterval = 1000.

#setup figure
plt.figure(2)
ax1 = plt.subplot(211)
ax2 = plt.subplot(212)

#Main Model Loop
while Time > EndTime:
    
    #Calculate integrated erosion for a tidal cycle
    Ez = np.zeros(len(Z))
    
    #bin the tidal elevations
    TidalWeights, TidalWaterLevels = np.histogram(SeaLevel+WaterLevels,bins=Z,density=True)

    for i in range(0,len(TidalWaterLevels)-1):
        
        if (TidalWaterLevels[i] < SeaLevel-TidalAmplitude): continue
        elif (TidalWaterLevels[i] > SeaLevel+TidalAmplitude): continue
            
        #Normalise elevations to water level for easy math
        Znormal = Z-(SeaLevel+TidalWaterLevels[i])
        
        #Find Position of Breaking Wave
        BreakingWaveInd = np.abs(Znormal+hb).argmin()          #Find bottom of surf zone
        
        if (X[BreakingWaveInd] == X[BreakingWaveInd-1]):
            XBreakingWave = X[BreakingWaveInd]
        elif (Znormal+hb)[BreakingWaveInd] < 0:
            XBreakingWave = X[BreakingWaveInd]+(X[BreakingWaveInd+1]-X[BreakingWaveInd])*((-hb-Znormal[BreakingWaveInd])/(Znormal[BreakingWaveInd+1]-Znormal[BreakingWaveInd]))
        elif (Znormal+hb)[BreakingWaveInd] > 0:
            XBreakingWave = X[BreakingWaveInd]+(X[BreakingWaveInd]-X[BreakingWaveInd-1])*((-hb-Znormal[BreakingWaveInd])/(Znormal[BreakingWaveInd]-Znormal[BreakingWaveInd-1]))
        else:
            print "problem"
            
        #Find Position of Shoreline
        ShorelineInd = np.abs(Znormal).argmin()                #Find shoreline
        
        if (X[ShorelineInd] == X[ShorelineInd-1]):
            XShoreline = X[ShorelineInd]
        elif Znormal[ShorelineInd] < 0:
            XShoreline = X[ShorelineInd]+(X[ShorelineInd+1]-X[ShorelineInd])*((0-Znormal[ShorelineInd])/(Znormal[ShorelineInd+1]-Znormal[ShorelineInd]))
        elif Znormal[ShorelineInd] > 0:
            XShoreline = X[ShorelineInd]+(X[ShorelineInd]-X[ShorelineInd-1])*((0-Znormal[ShorelineInd])/(Znormal[ShorelineInd]-Znormal[ShorelineInd-1]))
        else:
            XShoreline = X[ShorelineInd]
            
        #determine surf zone width
        SurfZoneWidth = XBreakingWave-XShoreline   #Width of surf zone
        
        #Calculate the surf force (kg/m^2) reaching the water line
        SurfForce = 0.5*rho_w*(Hb/0.78)*np.exp(-k*SurfZoneWidth) # - MinSurfForce
        
        #Calculate erosion across the platform
        ErosionMask = Z <= TidalWaterLevels[i]
        WaterDepth = TidalWaterLevels[i]-Z[ErosionMask]
        #Do erosion here including exponential decline with water depth
        Ez[ErosionMask] += TidalWeights[i]*M*SurfForce*np.exp(-WaterDepth)
    
    #Erode the coast, check for cliff overhangs and remove overhanging material
    XCliff = X[0]    
    for i in range(0,NoNodes):
        
        #Calculate
        #Do submarine erosion here
        X[i] = X[i] - Ez[i]
        
        #Assume cliff maintains vertical or sea ward sloping profile
        if (X[i] > XCliff):
            X[i] = XCliff
        
        XCliff = X[i]
            
    #PLOT?    
    if Time <= PlotTime:
        #do plotting here
        print Time
        PlotTime -= PlotInterval
        Color = Time/MaxTime
        ax1.plot(X,Z,'-',color=cm.Paired(Color))
        ax2.plot(Ez,Z,'-',color=cm.Paired(Color))
        
    #UPDATE TIME
    Time -= dt
    SeaLevel += SLR*dt
    
plt.show()