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

plt.figure(1,figsize=(3.75,2.2))

#plot data
plt.errorbar(Dist,N,yerr=Nerr,fmt='ko',mfc='white',zorder=10)
plt.xlabel ('Distance from cliff (m)')
plt.ylabel(r'\textsuperscript{10}Be Concentration (atoms g\textsuperscript{-1})')
plt.yticks([0,2500,5000,7500,10000,12500,15000])
plt.xlim(-10,350)
plt.ylim(0,15000)

plt.tight_layout()
plt.savefig("HG_CRN_results.pdf")

#plot model results
lc = [1.,0.67,0.28]
fc = [1.0,0.9,0.8]
Filename = "CRN_Model_Results_R10.0440_R20.0000_t0.0000_Beta0.0167_T2.4.txt"
X1, NModel1 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel1 += N[11]
plt.plot(X1,NModel1,'--',color=lc)

#plot model results
Filename = "CRN_Model_Results_R10.0480_R20.0000_t0.0000_Beta0.0167_T2.4.txt"
X2, NModel2 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel2 += N[11]
plt.plot(X2,NModel2,'-',color=lc)

#plot model results
Filename = "CRN_Model_Results_R10.0510_R20.0000_t0.0000_Beta0.0167_T2.4.txt"
X3, NModel3 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel3 += N[11]
plt.plot(X3,NModel3,'--',color=lc)

#fill in the curves
plt.fill(np.append(X1,X3[::-1]),np.append(NModel1,NModel3[::-1]),color=fc)
plt.plot([0,350],[N[11],N[11]],'k--')
plt.plot([0,0],[0,15000],'k-')
plt.plot([44,44],[0,15000],'k--')
plt.fill([0,350,350,0],[0,0,N[11],N[11]],color=[0.9,0.9,0.9])
plt.fill([0,44,44,0],[0,0,15000,15000],color=[0.9,0.9,0.9])
plt.savefig("HG_CRN_results_Single.pdf")
plt.show()

