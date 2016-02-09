# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:10:34 2016

@author: mhurst
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt

#Setup model domain
X = np.arange(-1000.,1001.,1.)
Z = np.zeros(len(X))
CliffPositionX = 1000.
CliffIndex = len(X)-1

MeanPlatformGradient = 1./100.
StepHeight = 0.5
RetreatRate = 0.2

#plot control
PlotTime = 0
PlotInterval = 1000
plt.figure(1)

#Loop through time
Time = 0
EndTime = 7000.
dt = 1.
while (Time < EndTime):
    #update time
    Time += dt
    
    #update cliff position
    CliffPositionX -= RetreatRate*dt
    if CliffPositionX < X[CliffIndex-1]:
        CliffIndex -= 1
    
    #Update the stepped profile
    for i in range(CliffIndex,len(X)):
        Z[i] = -(X[i]-CliffPositionX)*MeanPlatformGradient
    
    
    #Round to nearest step height
    Z = np.around(Z/StepHeight)*StepHeight

    if (Time >= PlotTime):
        PlotTime += PlotInterval
        plt.plot(CliffPositionX,0,'ro')
        plt.plot(X,Z,'k-')
        
plt.show()