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
HiroFile = HiroFolder + "Tide1m-Wave1m-Rock0.5-Weasd0.005-save_prof.txt"

#intiialise figure
plt.figure(1,figsize=(6,3))

#read file
f = open(HiroFile)
Lines = f.readlines()
NoLines = len(Lines)
NoNodes = len(Lines[0].strip().split(","))
Z = np.arange(-NoNodes/2, NoNodes/2)*0.1

#Setup time control
PlotTime = 0
PlotInterval = 1000

#read and plot each line
for i in range(0,NoLines):
    X = np.array(Lines[i].strip().split(",")).astype(np.float)
    if i == PlotTime or i == NoLines-1:
        print i,
        plt.plot(X,Z,'k-')
        PlotTime += PlotInterval
     
plt.xlabel("Distance (m)")
plt.ylabel("Elevation (m)")
plt.xlim(-10,np.max(X)+10)
plt.tight_layout()
plt.show()
