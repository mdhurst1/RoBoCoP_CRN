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

def make_plot(FileName,ColourMap):
    
    #create blank figure
    fig1 = plt.figure(1,figsize=(6,6))

    #First plot the morphology through time
    # declare the file and the axis
    ProfileName = FileName+"_ShoreProfile.xz"
    f = open(ProfileName,'r')
    Lines = f.readlines()
    NoLines = len(Lines)
    EndTime = 6000.
    ax1 = plt.subplot(211)

    #now plot CRN concentration through time
    # declare the file and the axis
    ConcentrationName = FileName+"_CRNConcentrations.xn"
    f2 = open(ConcentrationName,'r')
    Lines2 = f2.readlines()
    NoLines = len(Lines)
    ax2 = plt.subplot(212)

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
            ax1.plot(X,Z,'-',lw=1.5,color=ColourMap((Time+1000)/(EndTime+1000)))
            PlotTime += PlotInterval
        
        else:
            CliffPositionX = X[0]
            
            if (Time == PlotTime):
                ax1.plot(X,Z,'-',lw=1.5,color=ColourMap((Time+1000)/(EndTime+1000)))
                ax2.plot(X2,N,'-',lw=1.5,color=ColourMap((Time+1000)/(EndTime+1000)))
                PlotTime += PlotInterval
                
            TimeOld = Time
            CliffPositionXOld = CliffPositionX
    
    # tweak the plot
    #ax1.set_xticklabels([])
    ax1.set_xlabel("Distance (m)")
    ax1.set_ylabel("Elevation (m)")
    ax2.set_xlabel('Distance from Cliff (m)')
    ax2.set_ylabel('Concentration (atoms g$^{-1}$)')
    #ax1.set_xlim(-1300,100)
    #ax2.set_xlim(0,1000)
    ax2.set_ylim(0,20000)
    #ax2.set_ylim(-1800,100)
    
    plt.show()

if __name__ == "__main__":
    TideRanges = [1,2,4,8]
    ColourMaps = [cm.Reds,cm.Oranges,cm.Greens,cm.Purples]
    for i in range(0,len(TideRanges)):
        FileName = "../driver_files/MTR-"+str(TideRanges[i])+"m"
        ColourMap = ColourMaps[i]
        make_plot(FileName,ColourMap)
        plt.savefig(FileName+".pdf")
        plt.savefig(FileName+".png",dpi=200)
        plt.clf()