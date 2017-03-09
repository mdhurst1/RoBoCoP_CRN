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
    plt.figure(1,figsize=(6,4))

    #First plot the morphology through time
    # declare the file and the axis
    ProfileName = FileName+"ShoreProfile.xz"
    f = open(ProfileName,'r')
    Lines = f.readlines()
    NoLines = len(Lines)
    EndTime = NoLines
    
    # Only plot every 1 000 years
    PlotTime = 0
    PlotInterval = 1000
    
    ax1 = plt.subplot(111)

    #Get header info and setup X coord
    Header = Lines[0].strip().split(" ")
    for j in range(1,NoLines-1,2):
        
        Line = (Lines[j].strip().split(" "))
        Time = float(Line[0])
        
        #Read morphology
        X = np.array(Line[1:],dtype="float64")
        Z = np.array((Lines[j+1].strip().split(" "))[1:],dtype="float64")
        
        if (Time == PlotTime):
            ax1.plot(X,Z,'-',lw=1.5,color=ColourMap((Time+1000)/(EndTime+1000)))
            PlotTime += PlotInterval
                
    
    # tweak the plot
    #ax1.set_xticklabels([])
    ax1.set_xlabel("Distance (m)")
    ax1.set_ylabel("Elevation (m)")
    
    plt.show()

if __name__ == "__main__":
    FileName = "../driver_files/"
    ColourMap = cm.Reds
    make_plot(FileName,ColourMap)
        