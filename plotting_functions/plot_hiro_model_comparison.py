# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:10:34 2016

Script to plot the results of RoBoCoP and RockyCoastCRN

Martin Hurst,
March 7th 2016

@author: mhurst
"""

#import modules
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rc

# Customise figure style #
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
rc('font',size=8)
rc('ytick.major',pad=5)
rc('xtick.major',pad=5)
#rc('text', usetex=True)
padding = 5

def make_plots(HiroFileName, FileName, Parameters):
    
    #create blank figure
    plt.figure(1,figsize=(6.6,3.3))
    ax1 = plt.subplot(111)

    # Only plot every 1 000 years
    PlotTime = 0
    PlotInterval = 100
    EndTime = 1001
    
    #First plot the morphology through time
    # declare the file and the axis
    f = open(HiroFileName,'r')
    Lines = f.readlines()
    NoLines = len(Lines)
    f.close()
    
    #Get header info and setup X coord
    for j in range(0,NoLines):
        
        Line = (Lines[j].strip().split(","))
        Time = j
        #Read morphology
        X = 0.1*(np.array(Line[1:],dtype="float64") -1)
        NValues = len(X)
        Z = np.arange(-37.3,37.6, 0.1)
        print len(X), len(Z)
        if ((Time >= PlotTime) and (Time < EndTime)):
            ax1.plot(X,Z,'-',lw=1.5,color='r')
            #ax1.text(X[0],Z[0]-1,str(np.int(Time))+" years",rotation=-90)
            PlotTime += PlotInterval
            
    #First plot the morphology through time
    # declare the file and the axis
    f = open(FileName,'r')
    Lines = f.readlines()
    NoLines = len(Lines)
    f.close()
    
    # Get info on vertical from header
    Header = np.array(Lines[0].strip().split(" "),dtype=np.float)
    CliffHeight = Header[0]
    
    #Get header info and setup X coord
    PlotTime = 0
    for j in range(1,NoLines-1):
        
        Line = (Lines[j].strip().split(" "))
        Time = float(Line[0])
        
        #Read morphology
        X = np.array(Line[1:],dtype="float64")
        NValues = len(X)
        Z = np.linspace(CliffHeight,-CliffHeight, NValues)
        
        if ((Time >= PlotTime) and (Time < EndTime)):
            ax1.plot(X,Z,'-',lw=1.5,color=cm.coolwarm(Time/EndTime))
            ax1.text(X[0],Z[0]-1,str(np.int(Time))+" years",rotation=-90)
            PlotTime += PlotInterval
    
    plt.text(X[-1]+5,Z[-1]+1,(  "Tr = "+str(Parameters[1])+" m\n"+
                    "H = "+str(Parameters[2])+" m\n"+
                    "W = "+str(Parameters[3])+"\n"+
                    "R = "+str(Parameters[4])+"\n"+
                    "Br = "+str(Parameters[5])+"\n"+
                    "Bo = "+str(Parameters[6])),fontsize=8)

    # tweak the plot
    plt.xlabel("Distance (m)")
    plt.ylabel("Elevation (m)")
    xmin, xmax = ax1.get_xlim()
    plt.xlim(xmin-10,xmax+10)
    plt.ylim(-CliffHeight,CliffHeight)
    plt.tight_layout()
    plt.savefig(FileName.rstrip(".xz")+"_compare.png")
    plt.show()
    #plt.clf()

if __name__ == "__main__":
    HiroPath = "../results/Hiro_Test_Output2/"
    Path = "../results/Hiro_Ensemble_Runs/"
    
    #set parameter values explored
    Gradients = [90]
    TidalRanges = [1,4,8]
    WaveHeights = [1,3]
    WeatheringRates = [0.0001,0.001,0.01]
    Resistances = [0.001,0.01,0.1]
    BreakingCoefficients = [1]
    BrokenCoefficients = [1]

    #Loop across parameter space
    for a in range(0,len(Gradients)):
        for b in range(0,len(TidalRanges)):
            for c in range(0,len(WaveHeights)):
                for d in range(0,len(WeatheringRates)):
                    for e in range(0,len(Resistances)):
                        for f in range(0,len(BreakingCoefficients)):
                            for g in range(0,len(BrokenCoefficients)):

                                if WeatheringRates[d] == 0.0001:
                                    WR = "wea1"
                                elif WeatheringRates[d] == 0.001:
                                    WR = "wea3"
                                elif WeatheringRates[d] == 0.01:
                                    WR = "wea5"

                                HiroFileName = HiroPath+"P90-T"+str(TidalRanges[b])+"m-"+WR+"-co5-wave"+str(WaveHeights[c])+"m-resi"+str(Resistances[e])+"-save_prof.txt"
                                
                                if Gradients[a] == 90:
                                    Gradient = 0
                                Parameters = [Gradient,TidalRanges[b],WaveHeights[c],WeatheringRates[d],Resistances[e],BreakingCoefficients[f],BrokenCoefficients[g]]
                                FileName    = Path + ("ShoreProfile_G" + str(Gradient)
                                            + "_T_" + str(TidalRanges[b])
                                            + "_H_" + str(WaveHeights[c])
                                            + "_W_" + str(WeatheringRates[d])
                                            + "_R_" + str(Resistances[e])
                                            + "_Br_" + str(BreakingCoefficients[f])
                                            + "_Bo_" + str(BrokenCoefficients[g])
                                            + ".xz")
                                            
                                make_plots(HiroFileName, FileName, Parameters)
                                sys.exit("Temp break")