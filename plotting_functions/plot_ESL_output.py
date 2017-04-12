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
import matplotlib.animation as animation

# Customise figure style #
#rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
rc('font',size=8)
rc('ytick.major',pad=5)
rc('xtick.major',pad=5)
padding = 5

#create blank figure
fig = plt.figure(1,figsize=(6.6,5))

#First plot the morphology through time
# declare the file and the axis
FileName = "../results/ESL/"
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

ax1 = plt.axes([0.12,0.8,0.8,0.15])
ax2 = plt.axes([0.12,0.1,0.8,0.6])

#First plot the milankovitch curve
SeaLevelFileName = FileName+"SeaLevel.z"
T, SL = np.loadtxt(SeaLevelFileName, unpack=True)
T /= 1000
line, = ax1.plot(T,SL,'b-')
point, = ax1.plot(T[0],SL[0],'bo')
ax1.set_xlabel("Time (ka)")
ax1.set_ylabel("Eustatic Sea\nLevel (m)")
ax1.set_xlim(-700,10)
ax1.set_ylim(-70,70)


#plot the last timestep (this will never show but sets the extent)
Line = Lines[-1].strip().split(" ")
X = np.array(Line[1:],dtype="float64")
NValues = len(X)
Z = np.arange(CliffHeight,-CliffHeight-dz, -dz)
line, = ax2.plot(X,Z,'k-',lw=1.5)
ax2.set_xlabel("Distance (m)")
ax2.set_ylabel("Elevation (m)")
ax2.set_aspect(2)
ax2.set_xlim(ax2.get_xlim())
ax2.set_ylim(ax2.get_ylim())

SkyFillX = np.concatenate((np.array([-10,1000,1000]),X, np.array([-10])))
SkyFillZ = np.concatenate((np.array([100,100,10]),Z, np.array([-100])))
SkyFill, = ax2.fill(SkyFillX,SkyFillZ,color=[0.9,0.95,1.])

SLTideInd = np.argmin(np.abs(SL[0]-Z))
WaterFillX = np.concatenate((np.array([-10]),X[SLTideInd:],np.array([-10])))
WaterFillZ = np.concatenate((np.array([Z[SLTideInd]]),Z[SLTideInd:],np.array([Z[-1]])))
WaterFill, = plt.fill(WaterFillX,WaterFillZ,color=[0.7,0.85,1.])

#concatenate Z with cliff top vector
Z = np.arange(CliffHeight+dz,-CliffHeight-dz, -dz)

#print time
Time = int(T[-1]*1000)
text = ax2.text(450,-70,"Time is "+str(int(Time))+" yr")

f = open("filelist.txt","w")

##Loop across the data
for i in range(0,len(T)):
##for i in range(-10,0,1):
#    
    
    #update the point
    point.set_data(T[i], SL[i])

    #get the profile data
    Line = (Lines[i+1].strip().split(" "))
    Time = int(float(Line[0]))
    X = np.array(Line[1:],dtype="float64")
    X = np.concatenate([[10000],X])
    
    print i, T[i], Time
    
    line.set_data(X,Z)
    
    text.set_text("Time is "+str(Time)+" yr")
    
    
    SkyFillX = np.concatenate((np.array([-10,1000,1000]),X, np.array([-10])))
    SkyFillZ = np.concatenate((np.array([100,100,10]),Z, np.array([-100])))
    SkyFill.remove()
    SkyFill, = plt.fill(SkyFillX,SkyFillZ,color=[0.9,0.95,1.])
    
    SLTideInd = np.argmin(np.abs(SL[i]-Z))
    WaterFillX = np.concatenate((np.array([-10]),X[SLTideInd:],np.array([-10])))
    WaterFillZ = np.concatenate((np.array([Z[SLTideInd]]),Z[SLTideInd:],np.array([Z[-1]])))
    WaterFill.remove()
    WaterFill, = plt.fill(WaterFillX,WaterFillZ,color=[0.7,0.85,1.])
    
    FileName = 'test'+str(Time)+'.png'
    f.write(FileName+"\n")
    plt.savefig(FileName,dpi=100)

f.close()
plt.show()
#
#if __name__ == "__main__":
#    FileName = "../results/ESL/"
#    make_animation(FileName)
        