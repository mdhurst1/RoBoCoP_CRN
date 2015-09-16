# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 09:42:30 2015

@author: mhurst
"""

#IMPORT MODULES
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

#setup file name for HG
Filename = "HG_single_chain_MAD.out"

#Load the chain file
j, RetreatRate1, RetreatRate2, ChangeTime, BeachWidth, NewLike, LastLike, NAccepted, NRejected  = np.loadtxt(Filename, unpack=True, skiprows=2, usecols=(0,1,2,3,4,5,6,7,8))

#Most likely parameters Peak 1
FindRR1 = 0.046
FindBW = 40.1

MinV = 10000.
ind = 0

for i in range(0,len(j)):
    V = np.abs(RetreatRate1[i]-FindRR1)/FindRR1
    V += np.abs(BeachWidth[i]-FindBW)/FindBW
    if V<MinV:
        MinV = V
        ind = i

print "HG single Peak 1:", RetreatRate1[ind], ", BW",  BeachWidth[ind], ", -LL", -np.log(NewLike[ind])

#Most likely parameters Peak 2
FindRR1 = 0.082
FindBW = 15.3

MinV = 10000.
ind = 0

for i in range(0,len(j)):
    V = np.abs(RetreatRate1[i]-FindRR1)/FindRR1
    V += np.abs(BeachWidth[i]-FindBW)/FindBW
    if V<MinV:
        MinV = V
        ind = i

print "HG single Peak 2: RR1", RetreatRate1[ind], ", BW", BeachWidth[ind], ", -LL", -np.log(NewLike[ind])

#setup file name for HG
Filename = "HG_linear_chain_MAD.out"

#Load the chain file
i, RetreatRate1, RetreatRate2, ChangeTime, BeachWidth, NewLike, LastLike, NAccepted, NRejected  = np.loadtxt(Filename, unpack=True, skiprows=2, usecols=(0,1,2,3,4,5,6,7,8))

#Most likely parameters
FindRR1 = 0.053
FindRR2 = 0.059
FindBW = 44.3

MinV = 10000.
ind = 0

for i in range(0,len(j)):
    V = np.abs(RetreatRate1[i]-FindRR1)/FindRR1
    V += np.abs(RetreatRate2[i]-FindRR2)/FindRR2
    V += np.abs(BeachWidth[i]-FindBW)/FindBW
    
    if V<MinV:
        MinV = V
        ind = i

print "HG linear:", RetreatRate1[ind], ", RR2", RetreatRate2[ind],  ", BW", BeachWidth[ind], ", -LL", -np.log(NewLike[ind])

#setup file name for HG
Filename = "HG_change_chain_MAD.out"

#Load the chain file
i, RetreatRate1, RetreatRate2, ChangeTime, BeachWidth, NewLike, LastLike, NAccepted, NRejected  = np.loadtxt(Filename, unpack=True, skiprows=2, usecols=(0,1,2,3,4,5,6,7,8))

#Most likely parameters
FindRR1 = 0.055
FindRR2 = 0.032
FindCT = 408.5
FindBW = 43.7

MinV = 10000.
ind = 0

for i in range(0,len(j)):
    V = np.abs(RetreatRate1[i]-FindRR1)/FindRR1
    V += np.abs(RetreatRate2[i]-FindRR2)/FindRR2
    V += np.abs(ChangeTime[i]-FindCT)/FindCT
    V += np.abs(BeachWidth[i]-FindBW)/FindBW
    
    if V<MinV:
        MinV = V
        ind = i

print "HG change Peak 1:", RetreatRate1[ind], ", RR2", RetreatRate2[ind], ", CT", ChangeTime[ind], ", BW",  BeachWidth[ind], ", -LL", -np.log(NewLike[ind])

#Most likely parameters
FindRR1 = 0.153
FindRR2 = 0.062
FindCT = 2271.5
FindBW = 32.6

MinV = 10000.
ind = 0

for i in range(0,len(j)):
    V = np.abs(RetreatRate1[i]-FindRR1)/FindRR1
    V += np.abs(RetreatRate2[i]-FindRR2)/FindRR2
    V += np.abs(ChangeTime[i]-FindCT)/FindCT
    V += np.abs(BeachWidth[i]-FindBW)/FindBW
    
    if V<MinV:
        MinV = V
        ind = i

print "HG change Peak 2:", RetreatRate1[ind], ", RR2", RetreatRate2[ind], ", CT", ChangeTime[ind], ", BW",  BeachWidth[ind], ", -LL", -np.log(NewLike[ind])



#setup file name for HG
Filename = "BH_single_chain_MAD.out"

#Load the chain file
i, RetreatRate1, RetreatRate2, ChangeTime, BeachWidth, NewLike, LastLike, NAccepted, NRejected  = np.loadtxt(Filename, unpack=True, skiprows=2, usecols=(0,1,2,3,4,5,6,7,8))

#Most likely parameters
FindRR1 = 0.053
FindBW = 50.3

MinV = 10000.
ind = 0

for i in range(0,len(j)):
    V = np.abs(RetreatRate1[i]-FindRR1)/FindRR1
    V += np.abs(BeachWidth[i]-FindBW)/FindBW
    
    if V<MinV:
        MinV = V
        ind = i

print "BH single:", RetreatRate1[ind], ", BW",  BeachWidth[ind], ", -LL", -np.log(NewLike[ind])

#setup file name for HG
Filename = "BH_linear_chain_MAD.out"

#Load the chain file
i, RetreatRate1, RetreatRate2, ChangeTime, BeachWidth, NewLike, LastLike, NAccepted, NRejected  = np.loadtxt(Filename, unpack=True, skiprows=2, usecols=(0,1,2,3,4,5,6,7,8))

#Most likely parameters
FindRR1 = 0.109
FindRR2 = 0.032
FindBW = 52.6

MinV = 10000.
ind = 0

for i in range(0,len(j)):
    V = np.abs(RetreatRate1[i]-FindRR1)/FindRR1
    V += np.abs(RetreatRate2[i]-FindRR2)/FindRR2
    V += np.abs(BeachWidth[i]-FindBW)/FindBW
    
    if V<MinV:
        MinV = V
        ind = i

print "BH linear:", RetreatRate1[ind], ", RR2", RetreatRate2[ind], ", BW",  BeachWidth[ind], ", -LL", -np.log(NewLike[ind])

#setup file name for HG
Filename = "BH_change_chain_MAD.out"

#Load the chain file
i, RetreatRate1, RetreatRate2, ChangeTime, BeachWidth, NewLike, LastLike, NAccepted, NRejected  = np.loadtxt(Filename, unpack=True, skiprows=2, usecols=(0,1,2,3,4,5,6,7,8))

#Most likely parameters
FindRR1 = 0.044
FindRR2 = 0.257
FindCT = 278.4
FindBW = 41.0

MinV = 10000.
ind = 0

for i in range(0,len(j)):
    V = np.abs(RetreatRate1[i]-FindRR1)/FindRR1
    V += np.abs(RetreatRate2[i]-FindRR2)/FindRR2
    V += np.abs(ChangeTime[i]-FindCT)/FindCT
    V += np.abs(BeachWidth[i]-FindBW)/FindBW
    
    if V<MinV:
        MinV = V
        ind = i

print "BH change:", RetreatRate1[ind], ", RR2", RetreatRate2[ind], ", CT", ChangeTime[ind], ", BW",  BeachWidth[ind], ", -LL", -np.log(NewLike[ind])
