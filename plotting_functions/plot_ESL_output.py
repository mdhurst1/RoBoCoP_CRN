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
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
rc('font',size=8)
rc('ytick.major',pad=5)
rc('xtick.major',pad=5)
padding = 5

def make_animation(FileName):
    
    #create blank figure
    plt.figure(1,figsize=(6.6,6.6))

    #First plot the morphology through time
    # declare the file and the axis
    ProfileName = FileName+"ShoreProfile.xz"
    f = open(ProfileName,'r')
    Lines = f.readlines()
    NoLines = len(Lines)
    StartTime = float(Lines[1].strip().split(" ")[0])
    EndTime = float(Lines[-1].strip().split(" ")[0])

    # Get info on vertical from header
    Header = np.array(Lines[0].strip().split(" "),dtype=np.float)
    CliffHeight = Header[0]
    dz = Header[1]
    
    # Only plot every 1 000 years
    PlotTime = StartTime
    PlotInterval = 100000
    
    #ax1 = plt.axes([0.1,0.8,0.8,0.15])
    #ax2 = plt.axes([0.1,0.1,0.8,0.65])
#    
#    #Get header info and setup X coord
#    for j in range(1,NoLines-1):
#        
#        Line = (Lines[j].strip().split(" "))
#        Time = float(Line[0])
#        
#        #Read morphology
#        X = np.array(Line[1:],dtype="float64")
#        NValues = len(X)
#        Z = np.linspace(CliffHeight,-CliffHeight, NValues)
#            
#        HighTideInd = np.argmin(np.abs(0.5-Z)) 
#        SLTideInd = np.argmin(np.abs(0-Z))
#        LowTideInd = np.argmin(np.abs(-0.5-Z))
#        
#        SkyFillX = np.concatenate((np.array([-10,1000,1000]),X, np.array([-10])))
#        SkyFillZ = np.concatenate((np.array([100,100,10]),Z, np.array([-100])))
#        plt.fill(SkyFillX,SkyFillZ,color=[0.9,0.95,1.])
#        
#        WaterFillX = np.concatenate((np.array([-10]),X[HighTideInd:],np.array([-10])))
#        WaterFillZ = np.concatenate((np.array([Z[HighTideInd]]),Z[HighTideInd:],np.array([Z[-1]])))
#        plt.fill(WaterFillX,WaterFillZ,color=[0.7,0.85,1.])
#        
#        plt.plot([-10,X[HighTideInd]],[0.5,0.5],'--',color=[0.4,0.6,0.8])
#        plt.plot([-10,X[SLTideInd]],[0.,0.],'-',color=[0.4,0.6,0.8])
#        plt.plot([-10,X[LowTideInd]],[-0.5,-0.5],'--',color=[0.4,0.6,0.8])
#        
#        plt.plot(X,Z,'k-')
#        plt.plot([X[0],1000],[Z[0],Z[0]],'k-')
#            
#        if (Time >= PlotTime):
#            colour = Time/StartTime
#            ax1.plot(X,Z,'-',lw=1.5,color=cm.RdBu(colour))
#            PlotTime += PlotInterval
#    
#    ax1.plot(X,Z,'k-',lw=1.5)
#            
#        
#    print Z[0], Z[-1]           
#    # tweak the plot
#    #ax1.set_xticklabels([])
#    plt.xlabel("Distance (m)")
#    plt.ylabel("Elevation (m)")
#    xmin, xmax = ax1.get_xlim()
#    #plt.xlim(xmin-10,xmax+10)
#    #plt.ylim(-CliffHeight,CliffHeight)
#    #plt.ylim(-30,30)
#    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    FileName = "../results/ESL/"
    make_animation(FileName)
        