# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:10:34 2016

@author: mhurst
"""

#import modules
import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

#setup figure
fig = plt.figure(1,figsize=(6,6))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

# choose colour map
ColourMap = cm.copper

# List of erosion rates
Rates = [1,5,9]

for Rate in Rates:

	# get filename
	FileName ="../driver_files/RetreatRate_"+str(Rate)+"_test.xzn"

	# open the file
	f = open(FileName,'r')

	# get the lines and find out how many lines
	Lines = f.readlines()
	NoLines = len(Lines)

	#Get header info and setup X coord
	Header = Lines[0].strip().split(" ")
	NXNodes = float(Header[0])
	PlatformWidth = float(Header[1])

	MaxTime = float(Lines[1].strip().split(" ")[0])
    
	#Get data  
	XLine = Lines[-4].strip().split(" ")  
	ZLine = Lines[-2].strip().split(" ")
	NLine = Lines[-1].strip().split(" ")
	Time = float(ZLine[0])

	X = np.array(XLine[1:],dtype="float64")
	Z = np.array(ZLine[1:],dtype="float64")
	N = np.array(NLine[1:],dtype="float64")
    
	#mask for NDVs
	mask = Z != -9999
	Zplot = Z[mask]
	Xplot = X[mask]
	Nplot = N[mask]
    
	Colour = (Rate-np.min(Rates))/(np.max(Rates)-np.min(Rates))
	ax1.plot(Xplot,Zplot,'-',c=ColourMap(Colour))
	ax2.plot(Xplot,Nplot,'-',c=ColourMap(Colour))

#File2 = "../driver_files/Bideford_CRN.data"
#X,CRN,Error=np.loadtxt(File2,unpack=True,skiprows=1,usecols=(1,2,3),delimiter=" ")

#ax2.errorbar(X,CRN,fmt='o',yerr=Error,c='k')
#ax2.scatter(X,CRN)
    
ax2.set_xlabel("Distance (m)")
ax2.set_ylabel(r"$^{10}$Be Concentration")
ax1.set_ylabel("Elevation (m)")

#limit x axis to 250m
plt.xlim(0,250) 

plt.savefig("test_output_CB_3.png", dpi=300)  
