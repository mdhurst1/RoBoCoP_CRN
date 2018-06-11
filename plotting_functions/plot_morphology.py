# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:10:34 2016

@author: mhurst
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt

#open morphology file and read
FileName = "../driver_files/PlatformProfile.txt"
f = open(FileName,'r')
Lines = f.readlines()
NoLines = len(Lines)
#Get header info and setup X coord
Header = Lines[0].strip().split(" ")
NXNodes = float(Header[0])
PlatformWidth = float(Header[1])
X = np.linspace(0,PlatformWidth+1,NXNodes)
    
#setup figure
plt.figure(1,figsize=(6,3))

#Loop through data and plot
for i in range(1,NoLines,2):
    #Get data    
    Line = Lines[i].strip().split(" ")
    Time = float(Line[0])
    Z = np.array(Line[1:],dtype="float64")
    
    #mask for NDVs
    mask = Z != -9999
    Zplot = Z[mask]
    Xplot = X[mask]
    
    plt.plot(Xplot,Zplot,'k-')

plt.show()    
