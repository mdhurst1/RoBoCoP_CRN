# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:10:34 2016

Script to plot compare the results of RoBoCoP and RockyCoastCRN for transient 
and equillibrium runs

Martin Hurst,
June 25th 2016

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

# two figures, one for results of RoBoCoP runs, one for comparison to 
# equillibrium assumption
 
fig12 = plt.figure(12,figsize=(6,8))
fig13 = plt.figure(13,figsize=(6,6))

#First plot the morphology through time
# declare the file and the axis
FileName = "../results/RoBoCoP/No_SLR/RoBoCoP_ShoreProfile_NoSLR.xz"
f = open(FileName,'r')
Lines = f.readlines()
NoLines = len(Lines)
EndTime = 10000.
ax1 = fig12.add_axes([0.12,0.72,0.85,0.26])

#now plot CRN concentration through time
# declare the file and the axis
FileName = "../results/RoBoCoP/No_SLR/RoBoCoP_CRNConcentrations_NoSLR.xn"
f2 = open(FileName,'r')
Lines2 = f2.readlines()
NoLines = len(Lines)

#setup axes for plotting concentrations in both plots
ax2 = fig12.add_axes([0.12,0.37,0.85,0.3])
ax5 = fig13.add_axes([0.12,0.55,0.85,0.4])

#Setup the final axis for plotting retreat rate through time
ax3 = fig12.add_axes([0.12,0.06,0.35,0.24])

#setup axos for plotting platform gradient through time
ax4 = fig12.add_axes([0.62,0.06,0.35,0.24])

#setup axis on fig 13 for comparing peak concentrations
ax6 = fig13.add_axes([0.12,0.1,0.35,0.35])
ax7 = fig13.add_axes([0.62,0.1,0.35,0.35])

#Place holder lists for time, retreat rate and platform gradient
Times = []
RetreatRate = []
PlatformGradient = []
PeakN1 = []
PeakX1 = []

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
            ax1.plot(X,Z,'-',lw=1.5,color=cm.gray_r((Time+1000)/(EndTime+1000)))
            ax2.plot(X2-X2[0],N,'-',lw=1.5,color=cm.gray_r((Time+1000)/(EndTime+1000)))
            ax5.plot(X2-X2[0],N,'-',lw=1.5,color=cm.gray_r((Time+1000)/(EndTime+1000)))
            
            PeakN1.append(np.max(N))
            PeakX1.append(X2[np.argmax(N)]-X2[0])
            
            PlotTime += PlotInterval
            
        TimeOld = Time
        CliffPositionXOld = CliffPositionX

ax3.semilogy(np.array(Times)/1000,np.abs(RetreatRate),'k-')
ax4.plot(np.array(Times)/1000,np.abs(PlatformGradient),'k-')

print Times
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
ax1.set_xlim(-1300,100)
ax2.set_xlim(0,1000)
ax1.set_ylim(-4,4)

ax1.text(20,3.2,'(a)')
ax2.text(20,23000,'(b)')
ax3.text(0.6,4,'(c)')
ax4.text(8.8,0.01,'(d)')

RetreatRates = [0.3045, 0.1452, 0.0950, 0.0705, 0.0560, 0.0465, 0.0397, 0.0346, 0.0307, 0.0276]

PeakN2 = []
PeakX2 = []

for i in range (0,len(RetreatRates)):
    FileName = "../results/RoBoCoP/No_SLR/RetreatRate1_" + str(RetreatRates[i]) + ".xzn"
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
    
    PeakN2.append(np.max(N))
    PeakX2.append(X[np.argmax(N)])
    #setup color
    Color = (float(i)+1)/(float(len(RetreatRates))+1)
    ax5.plot(X,N,'--',color=cm.gray_r(Color),lw=1.5)

ax5.set_xlim(0,1000)
ax5.set_xlabel("Distance (m)")
ax5.set_ylabel('Concentration (atoms g$^{-1}$)')

print len(PeakN1), len(PeakN2)
for i in range(0,len(PeakN1)):
    Color = (float(i)+1)/(float(len(PeakN1))+1)
    print Color
    ax6.plot(PeakN1[i], PeakN2[i], 'o', mfc=cm.gray_r(Color),zorder=10)
    ax7.plot(PeakX1[i], PeakX2[i], 'o', mfc=cm.gray_r(Color))

ax6.plot([0,50000],[0,50000],'k--')
ax6.set_xlim(0,25000)
ax6.set_ylim(0,25000)
ax6.set_xlabel('Transient Concentration (atoms g$^{-1}$)')
ax6.set_ylabel('Equilibrium Concentration (atoms g$^{-1}$)')

ax7.plot([0,50000],[0,50000],'k--')
ax7.set_xlim(200,450)
ax7.set_ylim(200,450)
ax7.set_xlabel('Transient Position\nof Peak Concentration (m)')
ax7.set_ylabel('Equilibrium Position\nof Peak Concentration (m)')


ax5.text(20,23000,'(a)')
ax6.text(1500,22000,'(b)')
ax7.text(220,425,'(c)')

plt.show()