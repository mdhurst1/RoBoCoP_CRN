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
Filename = "HG_CRN_data_all.csv"
ID,Dist,Ncorr,Ncorrerr,N,Nerr = np.loadtxt(Filename,skiprows=1,delimiter=',',unpack=True,dtype=[('s5','S5'),('f4',float),('f5',float),('f6',float),('f7',float),('f8',float)])

plt.figure(1,figsize=(7.5,6.5))

#subplot 1
ax = plt.axes([0.2,0.71,0.36,0.27])

#plot data
plt.errorbar(Dist,N,yerr=Nerr,fmt='ko',mfc='white',zorder=10)

#plot model results
lc = [1.,0.67,0.28]
fc = [1.0,0.9,0.8]
Filename = "CRN_Model_Results_R10.0430_R20.0000_t0.0000_Beta0.0167_T2.4.txt"
X1, NModel1 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel1 += N[11]
plt.plot(X1,NModel1,'--',color=lc)

#plot model results
Filename = "CRN_Model_Results_R10.0460_R20.0000_t0.0000_Beta0.0167_T2.4.txt"
X2, NModel2 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel2 += N[11]
plt.plot(X2,NModel2,'-',color=lc)

#plot model results
Filename = "CRN_Model_Results_R10.0500_R20.0000_t0.0000_Beta0.0167_T2.4.txt"
X3, NModel3 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel3 += N[11]
plt.plot(X3,NModel3,'--',color=lc)

#fill in the curves
plt.fill(np.append(X1,X3[::-1]),np.append(NModel1,NModel3[::-1]),color=fc)

ax.set_xticklabels([])
plt.ylabel(r'\textsuperscript{10}Be Concentration (atoms g\textsuperscript{-1})')
plt.yticks([0,2500,5000,7500,10000,12500,15000])
plt.xlim(-10,350)
plt.ylim(0,15000)

#subplot 2
ax = plt.axes([0.6,0.71,0.36,0.27])
#plot data
plt.errorbar(Dist,N,yerr=Nerr,fmt='ko',mfc='white',zorder=10)

#plot model results
Filename = "CRN_Model_Results_R10.0710_R20.0000_t0.0000_Beta0.0167_T2.4.txt"
X1, NModel1 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel1 += N[11]
plt.plot(X1,NModel1,'--',color=lc)

#plot model results
Filename = "CRN_Model_Results_R10.0820_R20.0000_t0.0000_Beta0.0167_T2.4.txt"
X2, NModel2 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel2 += N[11]
plt.plot(X2,NModel2,'-',color=lc)

#plot model results
Filename = "CRN_Model_Results_R10.1000_R20.0000_t0.0000_Beta0.0167_T2.4.txt"
X3, NModel3 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel3 += N[11]
plt.plot(X3,NModel3,'--',color=lc)

#fill in the curves
plt.fill(np.append(X1,X3[::-1]),np.append(NModel1,NModel3[::-1]),color=fc)

ax.set_xticklabels([])
ax.set_yticklabels([])
plt.yticks([0,2500,5000,7500,10000,12500,15000])
plt.xlim(-10,350)
plt.ylim(0,15000)

#suplot 3
plt.axes([0.2,0.42,0.36,0.27])

#plot data
plt.errorbar(Dist,N,yerr=Nerr,fmt='ko',mfc='white',zorder=10)

#plot model results
lc = [0.5,0.6,1.0]
fc = [0.6,0.8,1.0]

Filename = "CRN_Model_Results_R10.0440_R20.0110_t104.1900_Beta0.0167_T2.4.txt"
X1, NModel1 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel1 += N[11]
plt.plot(X1,NModel1,'--',color=lc)

#plot model results
Filename = "CRN_Model_Results_R10.0550_R20.0320_t408.5000_Beta0.0167_T2.4.txt"
X2, NModel2 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel2 += N[11]
plt.plot(X2,NModel2,'-',color=lc)

#plot model results
Filename = "CRN_Model_Results_R10.0760_R20.1000_t1191.5000_Beta0.0167_T2.4.txt"
X3, NModel3 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel3 += N[11]
plt.plot(X3,NModel3,'--',color=lc)

#fill in the curves
plt.fill(np.append(X1,X3[::-1]),np.append(NModel1,NModel3[::-1]),color=fc)

plt.xlabel('Distance from Cliff (m)')
plt.ylabel(r'\textsuperscript{10}Be Concentration (atoms g\textsuperscript{-1})')
plt.yticks([0,2500,5000,7500,10000,12500,15000])
plt.xlim(-10,350)
plt.ylim(0,15000)

#subplot 4
ax = plt.axes([0.6,0.42,0.36,0.27])

#plot data
plt.errorbar(Dist,N,yerr=Nerr,fmt='ko',mfc='white',zorder=10)

#plot model results
lc = [0.5,0.6,1.0]
fc = [0.6,0.8,1.0]

Filename = "CRN_Model_Results_R10.1090_R20.0520_t1928.1000_Beta0.0167_T2.4.txt"
X1, NModel1 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel1 += N[11]
plt.plot(X1,NModel1,'--',color=lc)

#plot model results
Filename = "CRN_Model_Results_R10.1530_R20.0620_t2271.5000_Beta0.0167_T2.4.txt"
X2, NModel2 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel2 += N[11]
plt.plot(X2,NModel2,'-',color=lc)

#plot model results
Filename = "CRN_Model_Results_R10.1820_R20.0740_t2537.8000_Beta0.0167_T2.4.txt"
X3, NModel3 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel3 += N[11]
plt.plot(X3,NModel3,'--',color=lc)

#fill in the curves
plt.fill(np.append(X1,X3[::-1]),np.append(NModel1,NModel3[::-1]),color=fc)

plt.xlabel('Distance from Cliff (m)')
ax.set_yticklabels([])
plt.yticks([0,2500,5000,7500,10000,12500,15000])
plt.xlim(-10,350)
plt.ylim(0,15000)

#sunplot 5
plt.axes([0.4,0.08,0.36,0.27])
#plot data
plt.errorbar(Dist,N,yerr=Nerr,fmt='ko',mfc='white',zorder=10)
#plot model results
lc = [0.5,1.0,0.5]
fc = [0.8,1.0,0.8]

Filename = "CRN_Model_Results_R10.0500_R20.0520_t0.0000_Beta0.0167_T2.4.txt"
X1, NModel1 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel1 += N[11]
plt.plot(X1,NModel1,'--',color=lc)

#plot model results
Filename = "CRN_Model_Results_R10.0530_R20.0590_t0.0000_Beta0.0167_T2.4.txt"
X2, NModel2 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel2 += N[11]
plt.plot(X2,NModel2,'-',color=lc)

#plot model results
Filename = "CRN_Model_Results_R10.0600_R20.0680_t0.0000_Beta0.0167_T2.4.txt"
X3, NModel3 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel3 += N[11]
plt.plot(X3,NModel3,'--',color=lc)

#fill in the curves
plt.fill(np.append(X1,X3[::-1]),np.append(NModel1,NModel3[::-1]),color=fc)

plt.xlabel('Distance from Cliff (m)')
plt.ylabel(r'\textsuperscript{10}Be Concentration (atoms g\textsuperscript{-1})')
plt.yticks([0,2500,5000,7500,10000,12500,15000])
plt.xlim(-10,350)
plt.ylim(0,15000)

plt.savefig("HG_results_all.pdf")
plt.savefig("HG_results_all.png")
plt.show()
