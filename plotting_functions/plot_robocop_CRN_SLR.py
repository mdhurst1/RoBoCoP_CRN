# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:10:34 2016

Script to plot the results of RoBoCoP for the case of a transient model
run with relative sea level rise

Martin Hurst
University of Glasgow
July 18th 2016

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

#two figures, one for results, one for retreat rate and platform gradient
fig1 = plt.figure(1,figsize=(6,8))

#First plot the morphology through time
# declare the file and the axis
FileName = "../results/RoBoCoP/SLR/RoBoCoP_ShoreProfile_SLR.xz"
f = open(FileName,'r')
Lines = f.readlines()
NoLines = len(Lines)
EndTime = 10000.
ax1 = plt.axes([0.12,0.72,0.85,0.26])

#now plot CRN concentration through time
# declare the file and the axis
FileName = "../results/RoBoCoP/SLR/RoBoCoP_CRNConcentrations_SLR.xn"
f2 = open(FileName,'r')
Lines2 = f2.readlines()
NoLines = len(Lines)
ax2 = plt.axes([0.12,0.37,0.85,0.3])

#Setup the final axis for plotting retreat rate through time
ax3 = plt.axes([0.12,0.06,0.35,0.24])

#setup axos for plotting platform gradient through time
ax4 = plt.axes([0.62,0.06,0.35,0.24])

#Place holder lists for time, retreat rate and platform gradient
Times = []
RetreatRate = []
PlatformGradient = []

# Only plot every 1 000 years
PlotTime = 0
PlotInterval = 1000

#Get header info and setup X coord
Header = Lines[0].strip().split(" ")
for j in range(1,NoLines,2):
    
    Line = (Lines[j].strip().split(" "))
    Time = float(Line[0])
    
    #Read morphology
    X = np.array(Line[1:],dtype="float64")
    Z = np.array((Lines[j+1].strip().split(" "))[1:],dtype="float64")
    
    #Read concentrations
    X2 = np.array((Lines2[j].strip().split(" "))[1:],dtype="float64")
    N = np.array((Lines2[j+1].strip().split(" "))[1:],dtype="float64")
    
    if (j == 1):
        CliffPositionXOld = X[0]
        TimeOld = 0
        ax1.plot(X,Z,'-',lw=1.5,color=cm.gray_r((Time+1000)/(EndTime+1000)))
        PlotTime += PlotInterval
    
    else:
        CliffPositionX = X[0]
        Times.append(Time)
        RR = (CliffPositionXOld-CliffPositionX)/(TimeOld-Time)
        RetreatRate.append(RR)
        
        # Platform gradient, measure bewteen 0.5 and -0.9m elevation
        TopInd = np.argmin(np.abs(Z-0.5))
        BottomInd = np.argmin(np.abs(Z-(-0.9)))
        Grad = np.abs((Z[TopInd]-Z[BottomInd])/(X[TopInd]-X[BottomInd]))
        PlatformGradient.append(Grad)
        
        if (Time == PlotTime):
            ax1.plot(X,Z,'-',lw=1.5,color=cm.Blues((Time+1000)/(EndTime+1000)),label=int(Time/1000))
            ax2.plot(X2-X2[0],N,'-',lw=1.5,color=cm.Blues((Time+1000)/(EndTime+1000)))
            print Time, RR, np.median(RetreatRate), Grad
            PlotTime += PlotInterval
            
        TimeOld = Time
        CliffPositionXOld = CliffPositionX

ax3.loglog(np.array(Times)/1000,np.abs(RetreatRate),'k-')
ax4.plot(np.array(Times)/1000,np.abs(PlatformGradient),'k-')

# tweak the plot
#ax1.set_xticklabels([])
ax1.set_xlabel("Distance (m)")
ax1.set_ylabel("Elevation (m)")
ax2.set_xlabel('Distance from Cliff (m)')
ax2.set_ylabel('Concentration (atoms g$^{-1}$)')
ax3.set_xlabel('Time (ka)')
ax3.set_ylabel('RetreatRate (m yr$^{-1}$)')
ax4.set_xlabel('Time (ka)')
ax4.set_ylabel('Platform Gradient (m m$^{-1}$)')
ax1.set_xlim(-4100,100)
ax2.set_xlim(0,1000)
ax1.set_ylim(-5,10)
#ax2.set_ylim(-1800,100)

#Display legend
plt.rcParams.update({'legend.labelspacing':0.1}) 
plt.rcParams.update({'legend.columnspacing':1.0}) 
plt.rcParams.update({'legend.numpoints':1}) 
plt.rcParams.update({'legend.frameon':False}) 
plt.rcParams.update({'legend.handlelength':1.0}) 
plt.rcParams.update({'legend.fontsize':8})
ax1.legend(loc=3,ncol=2,title="time (K yr)")

#Load equilibrium results and add to ax2
FileName = "../results/RoBoCoP/SLR/Equilibrium_R0.32_SLR0.5.xzn"
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
    
#add to plot
ax2.plot(X,N,'k--',lw=1.5,label="Steady state model")
ax2.legend(loc=2,ncol=2)

plt.savefig("Transient_SLR.pdf")
plt.show()