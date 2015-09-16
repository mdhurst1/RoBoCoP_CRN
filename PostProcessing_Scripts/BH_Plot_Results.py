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

plt.figure(1,figsize=(3.75,2.2))
#plot data
plt.errorbar(Dist,N,yerr=Nerr,fmt='ko',mfc='white',zorder=12)
plt.ylabel(r'\textsuperscript{10}Be Concentration (atoms g\textsuperscript{-1})')
plt.yticks([0,2500,5000,7500,10000,12500,15000])
plt.xlim(0,300)
plt.ylim(0,15000)
plt.xlabel('Distance from cliff (m)')
plt.ylabel(r'\textsuperscript{10}Be Concnetration (atoms g\textsuperscript{-1}')
plt.tight_layout()

plt.savefig('BH_CRN_Results.pdf')

#plot model results
lc = [0.5,0.5,1.]
fc = [0.8,0.8,1.]
Filename = "CRN_Model_Results_R10.0140_R20.3000_t331.1000_Beta0.0167_T2.4.txt"
X1, NModel1 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel1 += N[1]
plt.plot(X1,NModel1,'--',color=lc,zorder=11)

#plot model results
Filename = "CRN_Model_Results_R10.0170_R20.2730_t430.4000_Beta0.0167_T2.4.txt"
X2, NModel2 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel2 += N[1]
plt.plot(X2,NModel2,'-',color=lc,zorder=11)

#plot model results
Filename = "CRN_Model_Results_R10.0260_R20.2000_t581.4500_Beta0.0167_T2.4.txt"
X3, NModel3 = np.loadtxt(Filename,skiprows=1,unpack=True)
NModel3 += N[1]
plt.plot(X3,NModel3,'--',color=lc,zorder=11)


#fill in the curves
plt.fill(np.append(X1,X3[::-1]),np.append(NModel1,NModel3[::-1]),color=fc,zorder=10)
plt.plot([0,300],[N[1],N[1]],'k--')
plt.plot([10,10],[0,15000],'k--')
plt.fill([0,10,10,0],[0,0,15000,15000],color=[0.9,0.9,0.9])
plt.fill([0,300,300,0],[0,0,N[1],N[1]],color=[0.9,0.9,0.9])

plt.savefig('BH_CRN_Results_Change.pdf')

plt.show()
