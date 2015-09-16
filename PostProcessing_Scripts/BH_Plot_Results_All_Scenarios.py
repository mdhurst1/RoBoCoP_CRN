# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 09:42:30 2015

@author: mhurst
"""

#IMPORT MODULES
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

##########################
# Customise figure style #
##########################
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
rc('font',size=8)
rc('ytick.major',pad=5)
rc('xtick.major',pad=5)
#rc('xtick', direction='out')
#rc('ytick', direction='out')
rc('text', usetex=True)
#rc('pdf',fonttype=42)
padding = 5

### load RAW data ###
Filename = "BH_CRN_data_all.csv"
ID,Dist,Ncorr,Ncorrerr,N,Nerr = np.loadtxt(Filename,skiprows=1,delimiter=',',unpack=True,dtype=[('s5','S5'),('f4',float),('f5',float),('f6',float),('f7',float),('f8',float)])

plt.figure(1,figsize=(3.75,5.5))
#subplot1
ax = plt.axes([0.25,0.7,0.7,0.27])
#plot data
plt.errorbar(Dist,N,yerr=Nerr,fmt='ko',mfc='white',zorder=10)

#plot model results
lc = [1.,0.67,0.28]
fc = [1.0,0.9,0.8]
Filename = "CRN_Model_Results_R10.0390_R20.0000_t0.0000_Beta0.0167_T2.4.txt"
X1, NModel1 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel1 += N[1]
plt.plot(X1,NModel1,'--',color=lc)

#plot model results
Filename = "CRN_Model_Results_R10.0530_R20.0000_t0.0000_Beta0.0167_T2.4.txt"
X2, NModel2 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel2 += N[1]
plt.plot(X2,NModel2,'-',color=lc)

#plot model results
Filename = "CRN_Model_Results_R10.0810_R20.0000_t0.0000_Beta0.0167_T2.4.txt"
X3, NModel3 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel3 += N[1]
plt.plot(X3,NModel3,'--',color=lc)

#fill in the curves
plt.fill(np.append(X1,X3[::-1]),np.append(NModel1,NModel3[::-1]),color=fc)

ax.set_xticklabels([])
plt.ylabel(r'\textsuperscript{10}Be Concentration (atoms g\textsuperscript{-1})')
#plt.yticks([4000,6000,8000,10000,12000])
plt.xlim(0,300)
plt.ylim(0,16000)


#subplot 2
ax = plt.axes([0.25,0.4,0.7,0.27])
#plot data
plt.errorbar(Dist,N,yerr=Nerr,fmt='ko',mfc='white',zorder=10)

#plot model results
lc = [0.5,0.6,1.0]
fc = [0.6,0.8,1.0]

Filename = "CRN_Model_Results_R10.0980_R20.3900_t108.8000_Beta0.0167_T2.4.txt"
X1, NModel1 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel1 += N[1]
plt.plot(X1,NModel1,'--',color=lc)

#plot model results
Filename = "CRN_Model_Results_R10.0440_R20.2570_t278.4000_Beta0.0167_T2.4.txt"
X2, NModel2 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel2 += N[1]
plt.plot(X2,NModel2,'-',color=lc)

#plot model results
Filename = "CRN_Model_Results_R10.0250_R20.1170_t657.6000_Beta0.0167_T2.4.txt"
X3, NModel3 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel3 += N[1]
plt.plot(X3,NModel3,'--',color=lc)

#fill in the curves
plt.fill(np.append(X1,X3[::-1]),np.append(NModel1,NModel3[::-1]),color=fc)

plt.ylabel(r'\textsuperscript{10}Be Concentration (atoms g\textsuperscript{-1})')
plt.xlim(0,300)
plt.ylim(0,16000)
ax.set_xticklabels([])

#subplot 3
plt.axes([0.25,0.1,0.7,0.27])
#plot data
plt.errorbar(Dist,N,yerr=Nerr,fmt='ko',mfc='white',zorder=10)

#plot model results
lc = [0.5,1.0,0.5]
fc = [0.8,1.0,0.8]

Filename = "CRN_Model_Results_R10.0440_R20.0110_t0.0000_Beta0.0167_T2.4.txt"
X1, NModel1 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel1 += N[1]
plt.plot(X1,NModel1,'--',color=lc)

#plot model results
Filename = "CRN_Model_Results_R10.1090_R20.0320_t0.0000_Beta0.0167_T2.4.txt"
X2, NModel2 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel2 += N[1]
plt.plot(X2,NModel2,'-',color=lc)

#plot model results
Filename = "CRN_Model_Results_R10.2140_R20.0960_t0.0000_Beta0.0167_T2.4.txt"
X3, NModel3 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel3 += N[1]
plt.plot(X3,NModel3,'--',color=lc)

plt.xlabel('Distance from Cliff (m)')
plt.ylabel(r'\textsuperscript{10}Be Concentration (atoms g\textsuperscript{-1})')
plt.xlim(0,300)
plt.ylim(0,16000)
plt.show()
