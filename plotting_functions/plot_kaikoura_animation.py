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
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rc

# Customise figure style #
rc('font',size=12)
rc('ytick.major',pad=5)
rc('xtick.major',pad=5)
#rc('font',**{'family':'serif','serif':['Helvetica']})
#rc('text', usetex=True)
padding = 5

#create blank figure
plt.figure(1,figsize=(16,9))

#First load the morphology through time
# declare the file and the axis
ProfileName = "../results/Kaikoura/Kaikoura_ShoreProfiles.xz"
f = open(ProfileName,'r')
MorphLines = f.readlines()
NoLines = len(MorphLines)
f.close()

dx = 0.1
dy = 0.1

EndTime = 10001

f = open("../results/Kaikoura/filelist.txt","w")

#Get figure extent
Line = MorphLines[5000].strip().split(" ")
#Read morphology
X = np.array(MorphLine[1:],dtype="float64")
print X

XMin = np.min(X)
XMax = np.max(X)
print XMin, XMax

#Loop through time and make each plot
j=1
for Time in range(0,EndTime+1,10):
    
    print Time, j
    
    #Get header info and setup X coord
    MorphLine = (MorphLines[j].strip().split(" "))
    j += 1
    
    #Read morphology
    X = np.array(MorphLine[1:],dtype="float64")
    Z = np.arange(10.,-10.01,-0.1)
        
    HighTideInd = np.argmin(np.abs(0.5-Z)) 
    SLTideInd = np.argmin(np.abs(0-Z))
    LowTideInd = np.argmin(np.abs(-0.5-Z))
    
    SkyFillX = np.concatenate((np.array([-10,1000,1000]),X, np.array([-10])))
    SkyFillZ = np.concatenate((np.array([100,100,10]),Z, np.array([-100])))
    plt.fill(SkyFillX,SkyFillZ,color=[0.9,0.95,1.])
    
    WaterFillX = np.concatenate((np.array([-10]),X[HighTideInd:],np.array([-10])))
    WaterFillZ = np.concatenate((np.array([Z[HighTideInd]]),Z[HighTideInd:],np.array([Z[-1]])))
    plt.fill(WaterFillX,WaterFillZ,color=[0.7,0.85,1.])
    
    plt.plot([-10,X[HighTideInd]],[0.5,0.5],'--',color=[0.4,0.6,0.8])
    plt.plot([-10,X[SLTideInd]],[0.,0.],'-',color=[0.4,0.6,0.8])
    plt.plot([-10,X[LowTideInd]],[-0.5,-0.5],'--',color=[0.4,0.6,0.8])
    
    plt.plot(X,Z,'k-')
    plt.plot([X[0],1000],[Z[0],Z[0]],'k-')
        
    # tweak the plot
    #ax1.set_xticklabels([])
    plt.xlabel(r"Distance (m)")
    plt.ylabel(r"Elevation (m)")
    plt.xlim(XMin,XMax)
    plt.ylim(-10,10)
    
    #write the time on the plot
    plt.text(-8,14,"Time in years: "+str(Time))
    
    #save plot and mark earthquake
    FigName = "../results/Kaikoura/Out"+str(Time)+".png"
    f.write(FigName+"\n")
    plt.savefig(FigName,dpi=100)
    if (Time == 5000):
        plt.text(10,9,'EARTHQUAKE')
        for i in range(0,26):
            FigName = "../results/Kaikoura/Out"+str(Time)+"_EQ"+str(i)+".png"
            f.write(FigName+"\n")
            plt.savefig(FigName,dpi=100)
            
    plt.clf()

f.close()

    
