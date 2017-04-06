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

def make_plot(Path, Parameters):
    
    #create blank figure
    plt.figure(1,figsize=(12,6))
    ax1 = plt.subplot(111)
    ax2 = plt.twinx()
    
    OffsetX = [0,500,800,0,600,200,350,900,0]
    OffsetZ = [50,70,90,0,30,40,80,30,120]
    
    for i in range(0,len(Parameters[0])):
        
        #Declare filenames
        MorphologyFileName = Path+"ShoreProfile_G1_T_"+str(Parameters[0][i])+"_H_"+str(Parameters[1][i])+"_W_"+str(Parameters[2][i])+"_R_"+str(Parameters[3][i])+"_Br_1_Bo_1.xz"
        ConcentrationsFileName = Path+"Concentrations_G1_T_"+str(Parameters[0][i])+"_H_"+str(Parameters[1][i])+"_W_"+str(Parameters[2][i])+"_R_"+str(Parameters[3][i])+"_Br_1_Bo_1.xn"
        
        #First plot the morphology through time
        # declare the file and the axis
        f = open(MorphologyFileName,'r')
        Lines = f.readlines()
        NoLines = len(Lines)
        EndTime = float(Lines[-1].strip().split(" ")[0])

        # Get info on vertical from header
        Header = np.array(Lines[0].strip().split(" "),dtype=np.float)
        CliffHeight = Header[0]
        dz = Header[1]
    
        Line = (Lines[-1].strip().split(" "))
        Time = float(Line[0])
        
        #Read morphology
        X = np.array(Line[1:],dtype="float64")+OffsetX[i]
        NValues = len(X)
        Z = np.arange(0,NValues)*-dz+CliffHeight+OffsetZ[i]
        
        #First plot the morphology through time
        # declare the file and the axis
        f = open(ConcentrationsFileName,'r')
        Lines = f.readlines()
        NoLines = len(Lines)
        dx = float(Lines[1].strip().split()[0])
        print dx
        #Read concentrations
        Line = (Lines[-1].strip().split(" "))
        N = np.array(Line[1:],dtype="float64")+OffsetZ[i]*1000
        Xn = np.arange(0,len(N))*dx + OffsetX[i]
        Mask = N != N[-1]
        N = N[Mask]
        Xn = Xn[Mask]
        
        ax1.plot(X,Z,'k-',lw=1.5)
        ax1.plot([X[0],X[-1]],[OffsetZ[i],OffsetZ[i]],'b-')
        ax2.plot(Xn,N,'r-',lw=1.5)
        
        ax1.text(X[-1]+20,Z[-1]-12,("Tr = "+str(Parameters[0][i])+" m\n"+
                    "H = "+str(Parameters[1][i])+" m\n"+
                    "W = "+str(Parameters[2][i])+"\n"+
                    "R = "+str(Parameters[3][i])),fontsize=8)

    # tweak the plot
    ax1.set_xlabel("Distance (m)")
    ax1.set_ylabel("Relative Elevation (m)")
    ax1.set_xlim(-100,1100)
    ax1.set_ylim(-25,150)
    ax2.set_ylim(-25000,150000)
    ax2.set_ylabel("Relative Concentration (atoms g$^{-1}$)")
    plt.tight_layout()
    plt.savefig("CRN_ensemble.png")
    plt.savefig("CRN_ensemble.pdf")
    #plt.clf()

if __name__ == "__main__":

    Path = "../results/Hiro_CRN_Ensembles/"
    
    #Loop across parameter space
    TidalRanges = [1,4,8,1,4,8,1,4,8]
    WaveHeights = [1,1,1,2,2,2,3,3,3]
    WeatheringRates = [0.01,0.001,0.01,0.0001,0.0001,0.001,0.0001,0.0001,0.001]
    Resistances = [0.1,0.1,0.1,0.001,0.001,0.01,0.01,0.1,0.001]
    Parameters = [TidalRanges,WaveHeights,WeatheringRates,Resistances]
    make_plot(Path,Parameters)