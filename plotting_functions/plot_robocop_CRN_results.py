# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:10:34 2016

Script to plot the results of coupled RoBoCoP and CRN simulations

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

FileName = "../data/RoBoCoP/ShoreProfile.xz"
f = open(FileName,'r')
ProfileLines = f.readlines()
NoLines = len(ProfileLines)

#now plot CRN concentration through time
FileName = "../data/RoBoCoP/CRNConcentrations.xn"
f = open(FileName,'r')
ConcentrationLines = f.readlines()
NoLines = len(ConcentrationLines)

#setup custom sized subplots
ax1 = plt.axes([0.12,0.675,0.86,0.3])
ax2 = plt.axes([0.12,0.1,0.86,0.55])

#Place holder lists for time and retreat rate
Times = []
RetreatRate = []
PrintTime = []
AverageRetreatRate = []
PlatformGradient = []

#Plotting control, plot every 1000 years
PlotTime = 1000.
PlotInterval = 1000.
EndTime = 5000.
TimeOld = 0

for j in range(1,NoLines,2):
    
    #Get time
    Line = (ProfileLines[j].strip().split(" "))
    Time = float(Line[0])
    
     #read in topographic profile
    X = np.array(Line[1:],dtype="float64")
    Z = np.array((ProfileLines[j+1].strip().split(" "))[1:],dtype="float64")
    
    #read in concentration data
    X2 = np.array((ConcentrationLines[j].strip().split(" "))[1:],dtype="float64")
    N = np.array((ConcentrationLines[j+1].strip().split(" "))[1:],dtype="float64")
    
    if j==1: 
        CliffPositionX = X[0]
        TimeOld = Time
        CliffPositionXOld = CliffPositionX
                   
        #plot profiles
        Color=0.2
        ax1.plot(X,Z,'-',color=cm.binary(Color), lw=1.5)
        datestring = '%.1f' % (Time/1000.)
        ax1.text(X[0]-4,6,datestring+" ka",rotation=270,color=cm.binary(Color))
        continue
    
    #only plot every plottime
    if (Time == PlotTime):
        #Get color for plotting
        Color = PlotTime/EndTime
        
        #plot profiles
        ax1.plot(X,Z,'-',color=cm.binary(Color), lw=1.5)
        datestring = '%.1f' % (Time/1000.)
        ax1.text(X[0]-4,6,datestring+" ka",rotation=270,color=cm.binary(Color))

        #plot concentrations
        ax2.plot(X2,N,'-',lw=1.5, color=cm.binary(Color))
        
        #update plot time
        PlotTime += PlotInterval        
    
        # Work out parameters to output
        # time
        PrintTime.append(Time)
        
        # Average retreat rate
        RR = (CliffPositionXOld-X[0])/(TimeOld-Time)
        AverageRetreatRate.append(RR)
        
        # Platform gradient, clip to cliff junction
        XClip = X[X>X[0]+1]
        ZClip = Z[X>X[0]+1]
        
        #Clip the last 200 m off too
        XClip = XClip[XClip<np.max(XClip)-200]
        ZClip = ZClip[XClip<np.max(XClip)-200]
        
        # find edge of platform as minimum gradient in X
        ind = np.argmin(np.gradient(XClip))
        ax1.plot(XClip[ind],ZClip[ind],'ro')
        
        #Calculate platform gradient
        PlatformGradient.append((ZClip[0]-ZClip[ind])/(XClip[0]-XClip[ind]))
    
    elif (Time > EndTime):
        break
    
    if j==1: 
        CliffPositionX = X[0]
        TimeOld = Time
        CliffPositionXOld = CliffPositionX
        continue
    
    CliffPositionX = X[0]
    Times.append(Time)
    RR = (CliffPositionXOld-CliffPositionX)/(TimeOld-Time)
    RetreatRate.append(RR)
    TimeOld = Time
    CliffPositionXOld = CliffPositionX
    
ax1.set_xlim(-1200,100)
ax2.set_xlim(-1200,100)
plt.xlabel('Distance (m)')
plt.ylabel('Elevation (m)')

ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.yaxis.set_ticks_position('left')
ax1.set_xticklabels([])
ax1.set_xticks([])

ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.yaxis.set_ticks_position('left')
ax2.xaxis.set_ticks_position('bottom')
#plt.subplot(212)
#plt.semilogy(Times,np.abs(RetreatRate),'k-')
#plt.xlabel("Time (years)")
#plt.ylabel("Retreat Rate (m/y)")
        
ax2.set_xlabel('Distance (m)')
ax2.set_ylabel('$^{10}$Be Concentration (atoms g$^{-1}$)')

#figure plotting retreat rates through time
ax3 = plt.axes([0.75,0.42,0.2,0.2])
ax3.loglog(Times,np.abs(RetreatRate),'k-')
ax3.set_xlabel('Time')
ax3.set_ylabel('RetreatRate (m yr$^{-1}$)')

ax1.text(-1170,8,"(a)")
ax2.text(-1170,7500,"(b)")
ax3.text(4000,3,"(c)")

#get average retreat rates at each time interval
MilleniaTimes = [1000,2000,3000,4000,5000]
AverageRate = [np.mean(RetreatRate[0:10]),np.mean(RetreatRate[10:20]),np.mean(RetreatRate[20:30]),np.mean(RetreatRate[30:40]),np.mean(RetreatRate[40:50])]

print PrintTime
print AverageRetreatRate
print PlatformGradient


#Plot equillibrium retreat results
RetreatRates = [0.146,0.089,0.76,0.073,0.072]
XNorms = [-700,-750,-750,-860,-870,-1000,-1100]

for i in range(0,len(RetreatRates)):
    FileName = "../driver_files/RetreatRate1_"+str(RetreatRates[i])+".xzn"
    f = open(FileName,'r')
    Lines = f.readlines()
    NoLines = len(Lines)
       
    X = np.array((Lines[NoLines-2].strip().split(" "))[1:],dtype="float64")
    N = np.array((Lines[NoLines-1].strip().split(" "))[1:],dtype="float64")
        
    Color = float(i)/float(len(RetreatRates))
    ax2.plot(X+XNorms[i],N,'-',color=cm.Reds(Color),label=str(RetreatRate[i]))
    

plt.savefig("../figures/RoBoCoP_CRN.pdf")
plt.show()