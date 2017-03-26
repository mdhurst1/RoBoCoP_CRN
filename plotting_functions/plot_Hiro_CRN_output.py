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
    plt.figure(1,figsize=(8,4))

    #First load the morphology through time
    # declare the file and the axis
    ProfileName = FileName+"ShoreProfile.xz"
    f = open(ProfileName,'r')
    MorphLines = f.readlines()
    NoLines = len(MorphLines)
    EndTime = float(MorphLines[-1].strip().split(" ")[0])
    f.close()
    
    #Second load CRN concentrations through time
    # declare the file and the axis
    ProfileName = FileName+"Concentrations.xn"
    f = open(ProfileName,'r')
    NLines = f.readlines()
    #NNoLines = len(NLines)
    #NEndTime = float(NLines[-1].strip().split(" ")[0])
    f.close()
    
    # Only plot every 1 000 years
    PlotTime = 0
    PlotInterval = 1000
    
    ax1 = plt.subplot(311)
    ax2 = plt.subplot(312)
    ax3 = plt.subplot(313)
    
    #Get header info and setup X coord
    for j in range(1,NoLines-1):
        
        MorphLine = (MorphLines[j].strip().split(" "))
        
        Time = float(N10Line[0])
        
        #Read morphology
        X = np.array(MorphLine[1:],dtype="float64")
        Z = np.arange(10.,-10.01,-0.1)
        X2 = np.arange(0,len(N10))*0.1
        
        if (Time == PlotTime):
            ax1.plot(X,Z,'-',lw=1.5,color=ColourMap((Time)/(EndTime)))
            PlotTime += PlotInterval
    
    N10Line = (NLines[-2].strip().split(" "))
    N14Line = (NLines[-1].strip().split(" "))
    N10 = np.array(N10Line[1:],dtype="float64")
    N14 = np.array(N14Line[1:],dtype="float64")
    mask = [N10!=N10[-1]]
    N10 = N10[mask]
    N14 = N14[mask]
    X2 = X2[mask]
           
    ax2.plot(X2,N10/1000.,'k-',lw=1.5)
    ax2.plot(X2,N14/1000.,'r-',lw=1.5)
    
    
    # tweak the plot
    #ax1.set_xticklabels([])
    
    ax1.set_ylabel("Elevation (m)")
    ax2.set_ylabel("Concentration (x 10$^3$ a g$^-1$)")
    xmin, xmax = ax1.get_xlim()
    ax1.set_xlim(xmin-10,xmax+10)
    ax2.set_xlim(xmin-10,xmax+10)
    ax1.set_ylim(-10,10)
    
    ax3.plot(X2,N10/N14)
    plt.xlabel("Distance (m)")
    ax3.set_ylabel("$^{10}$Be / $^{14}$C")
    
    plt.show()

if __name__ == "__main__":
    FileName = "../driver_files/"
    ColourMap = cm.gray_r
    make_plot(FileName,ColourMap)
        