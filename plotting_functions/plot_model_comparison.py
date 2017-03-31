# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 00:11:35 2017

Plot Hiro's model output and mine superimposed

@author: martin
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rc

# Customise figure style #
rc('font',size=12)
rc('ytick.major',pad=5)
rc('xtick.major',pad=5)
padding = 5

#load file
HiroFolder = "../results/Hiro_Test_Output/"
HiroFile = HiroFolder + "p90-1m-resi5-save_prof.txt"
CompareFile = HiroFolder + "ShoreProfile_T1_H1_W0.005.xz"
CompareFile = HiroFolder + "ShoreProfile.xz"

#intiialise figure
plt.figure(1,figsize=(6,3))

#read Hiro's model file
f = open(HiroFile)
Lines = f.readlines()
NoLines = len(Lines)
NoNodes = len(Lines[0].strip().split(","))
Z = np.arange(-NoNodes/2+1, NoNodes/2+1)*0.1

#Setup time control
PlotTime = 0
PlotInterval = 1000
EndTime = 6001

print "Hiro model",
#read and plot each line
for i in range(0,EndTime+1):
    
    if (i == PlotTime or i == NoLines-1) and (i<=NoLines-1):
        X = 0.1*(np.array(Lines[i].strip().split(",")).astype(np.float)-1)
        if X[0] == -1:
            break
        else:
            plt.plot(X,Z,'k-')
            plt.text(X[385]+0.1,Z[385],str(PlotTime))

        PlotTime += PlotInterval

plt.ylim(-10,10)

#Read my model file
f = open(CompareFile)
Lines = f.readlines()
NoLines = len(Lines)
NoNodes = len(Lines[1].strip().split(" "))-1
Z = np.arange(NoNodes/2, -NoNodes/2,-1)*0.1

#Setup time control
PlotTime = 0
PlotInterval = 1000

print "\nCompare",

#read and plot each line
for i in range(1,NoLines):
    Line = Lines[i].strip().split(" ")
    Time = int(Line[0])
    X = np.array(Line[1:]).astype(np.float)
    print Time, PlotTime
    if (Time == PlotTime) and (Time <= EndTime):
        print Time
        plt.plot(X,Z,'b-')
        PlotTime += PlotInterval
     
plt.xlabel("Distance (m)")
plt.ylabel("Elevation (m)")

plt.tight_layout()
plt.show()
