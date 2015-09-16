# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 09:42:30 2015

Script to plot the results of linear change in erosion rate for CRN coasts chain MCMC simulations

Martin Hurst, April 2014

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

##Set number of bins
Numbins = 100
CDF_Interpolate = [0.025,0.5,0.975]

#####################################
# PROCESS LINEAR CHAIN FOR HG FIRST #
#####################################
#
#Filename = "HG_linear_chain2.dout"
#print "Processing MCMC chain:", Filename
#
##Load the chain file
#i, RetreatRate1, RetreatRate2, ChangeTime, BeachWidth, ElevInit, NewLike, LastLike, NAccepted, NRejected  = np.loadtxt(Filename, unpack=True, skiprows=2, usecols=(0,1,2,3,4,5,6,7,8,9))
#
### Get accepted parameters for plotting convergence
#N = 0
#AcceptedMask = np.zeros(len(i),dtype=np.int)
#for j in range(0,len(i)):
#    if NAccepted[j] > N:
#        N += 1
#        AcceptedMask[j] = 1
#
#AcceptedRetreatRate1 = RetreatRate1[AcceptedMask == 1]
#AcceptedRetreatRate2 = RetreatRate2[AcceptedMask == 1]
#AcceptedChangeTime = ChangeTime[AcceptedMask == 1]
#AcceptedBeachWidth = BeachWidth[AcceptedMask == 1]
#Acceptedi = i[AcceptedMask == 1]
#

#
#### PROCESS AND PLOT EVERYTHING ###
#
###########################
## Customise figure style #
###########################
#from matplotlib import rc
##rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
#rc('font',size=8)
#rc('ytick.major',pad=5)
#rc('xtick.major',pad=5)
##rc('xtick', direction='out')
##rc('ytick', direction='out')
#rc('text', usetex=True)
##rc('pdf',fonttype=42)
#padding = 5
#
### function defining scientific notation
#def fmt(x, pos):
#    a, b = '{:.1e}'.format(x).split('e')
#    b = int(b)
#    return r'${} \times 10^{{{}}}$'.format(a, b)
#
##create dummy for blank axis labels
#blanklabels = []
#
## Now plot the progression of parameters through time
## plot the number of accepted and reject parameters through time too
#fig1 = plt.figure(1,figsize=(7.5,9.0))
#plt.subplots_adjust(0.1,0.1,0.9,0.9)
#ax = plt.subplot(411)
#plt.plot(i,RetreatRate1,'k-')
#plt.xlim(0,120000)
#ax.get_xaxis().set_ticks([])
#plt.ylabel('Retreat Rate 1 (m/y)')
#ax = plt.subplot(412)
#plt.plot(i,RetreatRate2,'k-')
#plt.xlim(0,120000)
#ax.get_xaxis().set_ticks([])
#plt.ylabel('Retreat Rate 2 (m/y)')
#ax = plt.subplot(413)
#plt.plot(i,BeachWidth,'k-')
#plt.xlim(0,120000)
#ax.get_xaxis().set_ticks([])
#plt.ylabel('Beach Width (m)')
#ax = plt.subplot(414)
#plt.plot(i,-np.log(NewLike),'k-')
#plt.xlim(0,120000)
#plt.xlabel('Number of Iterations')
#plt.ylabel('Negative Log-Likelihood')
#plt.suptitle('All Results')
#plt.tight_layout()
#plt.savefig("./HG_linear_chain_all.pdf")
#plt.savefig("./HG_linear_chain_all.png")
#
## Now plot the progression of parameters through time
## plot the number of accepted and reject parameters through time too
#fig2 = plt.figure(2,figsize=(7.5,9.))
#plt.subplots_adjust(0.1,0.1,0.9,0.9)
#
#ax = plt.subplot(411)
#plt.plot(Acceptedi,AcceptedRetreatRate1,'k-')
#ax.get_xaxis().set_ticks([])
#plt.ylim(0.0,0.25)
#plt.ylabel('Retreat Rate 1 (m/y)')
#plt.xlim(0,120000)
#
#plt.axes([0.15,0.81,0.2,0.05])
#plt.plot(Acceptedi[0:1000],AcceptedRetreatRate1[0:1000],'k-')
#plt.xlim(0,1000)
#plt.yticks([0.1,0.12,0.14,0.16,0.18,0.2])
#
#ax = plt.subplot(412)
#plt.plot(Acceptedi,AcceptedRetreatRate2,'k-')
#ax.xaxis.set_ticklabels([])
#plt.ylim(0.01,0.10)
#plt.ylabel('Retreat Rate 2 (m/y)')
#plt.xlim(0,120000)
#
#plt.axes([0.15,0.68,0.2,0.05])
#plt.plot(Acceptedi[0:1000],AcceptedRetreatRate2[0:1000],'k-')
#plt.xlim(0,1000)
#plt.yticks([0.02,0.04,0.06,0.08,0.1])
#
#ax = plt.subplot(413)
#plt.plot(Acceptedi,AcceptedBeachWidth,'k-')
#ax.xaxis.set_ticklabels([])
#plt.ylabel('Beach Width (m)')
#plt.xlim(0,120000)
#
#plt.axes([0.15,0.32,0.2,0.05])
#plt.plot(Acceptedi[0:1000],AcceptedBeachWidth[0:1000],'k-')
#plt.xlim(0,1000)
#plt.yticks([0,20,40,60])
#
#ax = plt.subplot(414)
#plt.plot(i,-np.log(LastLike),'k-')
#plt.ylim(40,60)
#plt.xlabel('Number of Iterations')
#plt.ylabel('Negative Log-Likelihood')
#plt.xlim(0,120000)
#
#plt.axes([0.15,0.18,0.2,0.05])
#plt.plot(i[0:1000],-np.log(LastLike[0:1000]),'k-')
#plt.xlim(0,1000)
#plt.yticks([40,50,60,70,80])
#
#plt.suptitle('Accepted Results')
#plt.tight_layout()
#plt.savefig("./HG_linear_chain_accepted.pdf")
#plt.savefig("./HG_linear_chain_accepted.png")
#
############ DO PDF AND CDF plots ############
#
##First for Retreat Rate 1
#Values, Edges = np.histogram(RetreatRate1,Numbins,normed=True,weights=NewLike)
#Values /= np.sum(Values)
#CDF = np.cumsum(Values)
#LeftEdges = Edges[:-1]
#RightEdges = Edges[1:]
#Midpoints = (RightEdges + LeftEdges)/2
#Result = np.interp(CDF_Interpolate,CDF,Midpoints)
#print "Retreat Rate 1:", Result
#LikelyRR1 = Result[1]
#
#fig3 = plt.figure(3,figsize=(7.5,9))
#plt.subplot(321)
#plt.bar(Edges[:-1],Values,width = Edges[1]-Edges[0])
#plt.xlim(0.1,0.25)
#plt.xlabel('Retreat Rate 1 (m/y)')
#plt.ylabel('Probability Density')
#
#plt.subplot(322)
#plt.plot(Midpoints,CDF,'k-')
#MinInd = np.argmin(np.abs(CDF-0.025))
#MaxInd = np.argmin(np.abs(CDF-0.975))
#MedInd = np.argmin(np.abs(CDF-0.5))
#plt.plot([0,Result[0]],[CDF_Interpolate[0],CDF_Interpolate[0]],'r--')
#plt.plot([0,Result[2]],[CDF_Interpolate[2],CDF_Interpolate[2]],'r--')
#plt.plot([0,Result[1]],[CDF_Interpolate[1],CDF_Interpolate[1]],'r-')
#
#plt.fill(np.append(Midpoints[MinInd:MaxInd+1],[Result[2],Result[2],Result[0]]), np.append(CDF[MinInd:MaxInd+1],[CDF_Interpolate[2],0,0]), color=[0.8,0.8,0.8])
#plt.plot([Result[2],Result[2]],[CDF_Interpolate[2],0],'k--')
#plt.plot([Result[0],Result[0]],[CDF_Interpolate[0],0],'k--')
#plt.plot([Result[1],Result[1]],[CDF_Interpolate[1],0],'k-')
#plt.xlim(0.1,0.25)
#plt.ylim(0,1)
#plt.xlabel('Retreat Rate 1 (m/y)')
#plt.ylabel('Cumulative Probability Density')
#
##Retreat Rate 2
#Values, Edges = np.histogram(RetreatRate2,Numbins,normed=True,weights=NewLike)
#Values /= np.sum(Values)
#CDF = np.cumsum(Values)
#LeftEdges = Edges[:-1]
#RightEdges = Edges[1:]
#Midpoints = (RightEdges + LeftEdges)/2
#Result = np.interp(CDF_Interpolate,CDF,Midpoints)
#print "Retreat Rate 2:", Result
#LikelyRR2 = Result[1]
#
#plt.subplot(323)
#plt.bar(Edges[:-1],Values,width = Edges[1]-Edges[0])
#plt.xlim(0.01,0.06)
#plt.xlabel('Retreat Rate 2 (m/y)')
#plt.ylabel('Probability Density')
#
#plt.subplot(324)
#plt.plot(Midpoints,CDF,'k-')
#MinInd = np.argmin(np.abs(CDF-0.025))
#MaxInd = np.argmin(np.abs(CDF-0.975))
#MedInd = np.argmin(np.abs(CDF-0.5))
#plt.plot([0,Result[0]],[CDF_Interpolate[0],CDF_Interpolate[0]],'r--')
#plt.plot([0,Result[2]],[CDF_Interpolate[2],CDF_Interpolate[2]],'r--')
#plt.plot([0,Result[1]],[CDF_Interpolate[1],CDF_Interpolate[1]],'r-')
#
#plt.fill(np.append(Midpoints[MinInd:MaxInd+1],[Result[2],Result[2],Result[0]]), np.append(CDF[MinInd:MaxInd+1],[CDF_Interpolate[2],0,0]), color=[0.8,0.8,0.8])
#plt.plot([Result[2],Result[2]],[CDF_Interpolate[2],0],'k--')
#plt.plot([Result[0],Result[0]],[CDF_Interpolate[0],0],'k--')
#plt.plot([Result[1],Result[1]],[CDF_Interpolate[1],0],'k-')
#plt.xlim(0.01,0.06)
#plt.ylim(0,1)
#plt.xlabel('Retreat Rate 2 (m/y)')
#plt.ylabel('Cumulative Probability Density')
#
#
##Finally for beach width
#Values, Edges = np.histogram(BeachWidth,Numbins,normed=True,weights=NewLike)
#Values /= np.sum(Values)
#CDF = np.cumsum(Values)
#LeftEdges = Edges[:-1]
#RightEdges = Edges[1:]
#Midpoints = (RightEdges + LeftEdges)/2
#Result = np.interp(CDF_Interpolate,CDF,Midpoints)
#print "Beach Width:", Result
#LikelyBW = Result[1]
#
#plt.subplot(325)
#plt.bar(Edges[:-1],Values,width = Edges[1]-Edges[0])
#plt.xlim(30,60)
#plt.xlabel('BeachWidth (m)')
#plt.ylabel('Probability Density')
#
#plt.subplot(326)
#plt.plot(Midpoints,CDF,'k-')
#MinInd = np.argmin(np.abs(CDF-0.025))
#MaxInd = np.argmin(np.abs(CDF-0.975))
#MedInd = np.argmin(np.abs(CDF-0.5))
#plt.plot([0,Result[0]],[CDF_Interpolate[0],CDF_Interpolate[0]],'r--')
#plt.plot([0,Result[2]],[CDF_Interpolate[2],CDF_Interpolate[2]],'r--')
#plt.plot([0,Result[1]],[CDF_Interpolate[1],CDF_Interpolate[1]],'r-')
#
#plt.fill(np.append(Midpoints[MinInd:MaxInd+1],[Result[2],Result[2],Result[0]]), np.append(CDF[MinInd:MaxInd+1],[CDF_Interpolate[2],0,0]), color=[0.8,0.8,0.8])
#plt.plot([Result[2],Result[2]],[CDF_Interpolate[2],0],'k--')
#plt.plot([Result[0],Result[0]],[CDF_Interpolate[0],0],'k--')
#plt.plot([Result[1],Result[1]],[CDF_Interpolate[1],0],'k-')
#plt.xlim(30,60)
#plt.ylim(0,1)
#plt.xlabel('Beach Width (m)')
#plt.ylabel('Cumulative Probability Density')
#
#plt.tight_layout()
#plt.savefig('./HG_linear_chain_PDFs.pdf')
#plt.savefig('./HG_linear_chain_PDFs.png')
#
#ind = np.argmax(NewLike)
#print "Most Likely RR1", RetreatRate1[ind], "RR2", RetreatRate2[ind], "BW", BeachWidth[ind]
#print "Likelihood", NewLike[ind]
#print "-Log Likelihood", -np.log(NewLike[ind])
#

#####################################################
#       NOW PLOT change CHAIN FOR BEACHY HEAD       #
#####################################################

#setup file name for BH

Filename = "BH_linear_chain2.dout"
print "Processing MCMC chain:", Filename

#Load the chain file
i, RetreatRate1, RetreatRate2, ChangeTime, BeachWidth, ElevInit, NewLike, LastLike, NAccepted, NRejected  = np.loadtxt(Filename, unpack=True, skiprows=2, usecols=(0,1,2,3,4,5,6,7,8,9))

# Get accepted parameters for plotting convergence
N = 0
AcceptedMask = np.zeros(len(i),dtype=np.int)
for j in range(0,len(i)):
    if NAccepted[j] > N:
        N += 1
        AcceptedMask[j] = 1

AcceptedRetreatRate1 = RetreatRate1[AcceptedMask == 1]
AcceptedRetreatRate2 = RetreatRate2[AcceptedMask == 1]
AcceptedChangeTime = ChangeTime[AcceptedMask == 1]
AcceptedBeachWidth = BeachWidth[AcceptedMask == 1]
Acceptedi = i[AcceptedMask == 1]

## PROCESS AND PLOT EVERYTHING ###

# Now plot the progression of parameters through time
# plot the number of accepted and reject parameters through time too
fig4 = plt.figure(4,figsize=(7.5,9.0))
plt.subplots_adjust(0.1,0.1,0.9,0.9)
ax = plt.subplot(511)
plt.plot(i,RetreatRate1,'k-')
ax.xaxis.set_ticklabels([])
plt.xlabel('Number of Iterations')
plt.ylabel('Retreat Rate 1 (m/y)')
ax = plt.subplot(512)
plt.plot(i,RetreatRate2,'k-')
ax.xaxis.set_ticklabels([])
plt.xlabel('Number of Iterations')
plt.ylabel('Retreat Rate 2 (m/y)')
ax = plt.subplot(513)
plt.plot(i,ChangeTime,'k-')
ax.xaxis.set_ticklabels([])
plt.xlabel('Number of Iterations')
plt.ylabel('Change Time (years BP)')
ax = plt.subplot(514)
plt.plot(i,BeachWidth,'k-')
ax.xaxis.set_ticklabels([])
plt.xlabel('Number of Iterations')
plt.ylabel('Beach Width (m)')
plt.subplot(515)
plt.plot(i,np.abs(np.log(NewLike)),'k-')
plt.xlabel('Number of Iterations')
plt.ylabel('Negative Log-Likelihood')
plt.suptitle('All Results')
plt.tight_layout()
plt.savefig("./BH_linear_chain_all.pdf")
plt.savefig("./BH_linear_chain_all.png")

# Now plot the progression of parameters through time
# plot the number of accepted and reject parameters through time too
fig5 = plt.figure(5,figsize=(7.5,9.))
plt.subplots_adjust(0.1,0.1,0.9,0.9)

ax = plt.subplot(411)
plt.plot(Acceptedi,AcceptedRetreatRate1,'k-')
ax.get_xaxis().set_ticks([])
#plt.ylim(0.04,0.12)
plt.ylabel('Retreat Rate 1 (m/y)')

plt.axes([0.15,0.9,0.22,0.05])
plt.plot(Acceptedi[0:500],AcceptedRetreatRate1[0:500],'k-')
plt.xlim(0,1000)
plt.yticks([0.02,0.04,0.06,0.08,0.1])

ax = plt.subplot(412)
plt.plot(Acceptedi,AcceptedRetreatRate2,'k-')
ax.xaxis.set_ticklabels([])
#plt.ylim(0.00,0.16)
plt.ylabel('Retreat Rate 2 (m/y)')

plt.axes([0.15,0.67,0.2,0.05])
plt.plot(Acceptedi[0:500],AcceptedRetreatRate2[0:500],'k-')
plt.xlim(0,1000)
plt.yticks([0.02,0.04,0.06,0.08,0.1])

ax = plt.subplot(413)
plt.plot(Acceptedi,AcceptedBeachWidth,'k-')
ax.xaxis.set_ticklabels([])
plt.ylim(0,60)
plt.ylabel('Beach Width (m)')

plt.axes([0.15,0.32,0.2,0.05])
plt.plot(Acceptedi[0:500],AcceptedBeachWidth[0:500],'k-')
plt.xlim(0,1000)
plt.yticks([20,40,60,80,100])

ax = plt.subplot(414)
plt.plot(i,np.abs(np.log(LastLike)),'k-')
plt.xlabel('Number of Iterations')
plt.ylabel('Negative Log-Likelihood')

plt.axes([0.15,0.18,0.2,0.05])
plt.plot(i[0:1000],np.abs(np.log(LastLike[0:1000])),'k-')
plt.xlim(0,1000)
plt.yticks([120,125,130,135,140])

plt.suptitle('Accepted Results')
plt.tight_layout()
plt.savefig("./BH_linear_chain_accepted.pdf")
plt.savefig("./BH_linear_chain_accepted.png")

########### DO PDF AND CDF plots ############

#First for Retreat Rate 1
Values, Edges = np.histogram(RetreatRate1,Numbins,normed=True,weights=NewLike)
Values /= np.sum(Values)
CDF = np.cumsum(Values)
LeftEdges = Edges[:-1]
RightEdges = Edges[1:]
Midpoints = (RightEdges + LeftEdges)/2
Result = np.interp(CDF_Interpolate,CDF,Midpoints)
print "Retreat Rate 1:", Result
LikelyRR1 = Result[1]

fig6 = plt.figure(6,figsize=(7.5,7.5))
plt.subplot(321)
plt.bar(Edges[:-1],Values,width = Edges[1]-Edges[0])
plt.xlim(0.01,0.05)
plt.xlabel('Retreat Rate 1 (m/y)')
plt.ylabel('Probability Density')

plt.subplot(322)
plt.plot(Midpoints,CDF,'k-')
MinInd = np.argmin(np.abs(CDF-0.025))
MaxInd = np.argmin(np.abs(CDF-0.975))
MedInd = np.argmin(np.abs(CDF-0.5))
plt.plot([0,Result[0]],[CDF_Interpolate[0],CDF_Interpolate[0]],'r--')
plt.plot([0,Result[2]],[CDF_Interpolate[2],CDF_Interpolate[2]],'r--')
plt.plot([0,Result[1]],[CDF_Interpolate[1],CDF_Interpolate[1]],'r-')

plt.fill(np.append(Midpoints[MinInd:MaxInd+1],[Result[2],Result[2],Result[0]]), np.append(CDF[MinInd:MaxInd+1],[CDF_Interpolate[2],0,0]), color=[0.8,0.8,0.8])
plt.plot([Result[2],Result[2]],[CDF_Interpolate[2],0],'k--')
plt.plot([Result[0],Result[0]],[CDF_Interpolate[0],0],'k--')
plt.plot([Result[1],Result[1]],[CDF_Interpolate[1],0],'k-')
plt.xlim(0.01,0.05)
plt.ylim(0,1)
plt.xlabel('Retreat Rate 1 (m/y)')
plt.ylabel('Cumulative Probability Density')

#Retreat Rate 2
Values, Edges = np.histogram(RetreatRate2,Numbins,normed=True,weights=NewLike)
Values /= np.sum(Values)
CDF = np.cumsum(Values)
LeftEdges = Edges[:-1]
RightEdges = Edges[1:]
Midpoints = (RightEdges + LeftEdges)/2
Result = np.interp(CDF_Interpolate,CDF,Midpoints)
print "Retreat Rate 2:", Result
LikelyRR2 = Result[1]

plt.subplot(323)
plt.bar(Edges[:-1],Values,width = Edges[1]-Edges[0])
plt.xlim(0.001,0.1)
plt.xlabel('Retreat Rate 2 (m/y)')
plt.ylabel('Probability Density')

plt.subplot(324)
plt.plot(Midpoints,CDF,'k-')
MinInd = np.argmin(np.abs(CDF-0.025))
MaxInd = np.argmin(np.abs(CDF-0.975))
MedInd = np.argmin(np.abs(CDF-0.5))
plt.plot([0,Result[0]],[CDF_Interpolate[0],CDF_Interpolate[0]],'r--')
plt.plot([0,Result[2]],[CDF_Interpolate[2],CDF_Interpolate[2]],'r--')
plt.plot([0,Result[1]],[CDF_Interpolate[1],CDF_Interpolate[1]],'r-')

plt.fill(np.append(Midpoints[MinInd:MaxInd+1],[Result[2],Result[2],Result[0]]), np.append(CDF[MinInd:MaxInd+1],[CDF_Interpolate[2],0,0]), color=[0.8,0.8,0.8])
plt.plot([Result[2],Result[2]],[CDF_Interpolate[2],0],'k--')
plt.plot([Result[0],Result[0]],[CDF_Interpolate[0],0],'k--')
plt.plot([Result[1],Result[1]],[CDF_Interpolate[1],0],'k-')
plt.xlim(0.001,0.1)
plt.ylim(0,1)
plt.xlabel('Retreat Rate 2 (m/y)')
plt.ylabel('Cumulative Probability Density')

#Now for Change Time
Values, Edges = np.histogram(ChangeTime,Numbins,normed=True,weights=NewLike)
Values /= np.sum(Values)
CDF = np.cumsum(Values)
LeftEdges = Edges[:-1]
RightEdges = Edges[1:]
Midpoints = (RightEdges + LeftEdges)/2
Result = np.interp(CDF_Interpolate,CDF,Midpoints)
print "Change Time:", Result

#Finally for beach width
Values, Edges = np.histogram(BeachWidth,Numbins,normed=True,weights=NewLike)
Values /= np.sum(Values)
CDF = np.cumsum(Values)
LeftEdges = Edges[:-1]
RightEdges = Edges[1:]
Midpoints = (RightEdges + LeftEdges)/2
Result = np.interp(CDF_Interpolate,CDF,Midpoints)
print "Beach Width:", Result
LikelyBW = Result[1]


plt.subplot(325)
plt.bar(Edges[:-1],Values,width = Edges[1]-Edges[0])
plt.xlim(20,100)
plt.xlabel('BeachWidth (m)')
plt.ylabel('Probability Density')

plt.subplot(326)
plt.plot(Midpoints,CDF,'k-')
MinInd = np.argmin(np.abs(CDF-0.025))
MaxInd = np.argmin(np.abs(CDF-0.975))
MedInd = np.argmin(np.abs(CDF-0.5))
plt.plot([0,Result[0]],[CDF_Interpolate[0],CDF_Interpolate[0]],'r--')
plt.plot([0,Result[2]],[CDF_Interpolate[2],CDF_Interpolate[2]],'r--')
plt.plot([0,Result[1]],[CDF_Interpolate[1],CDF_Interpolate[1]],'r-')

plt.fill(np.append(Midpoints[MinInd:MaxInd+1],[Result[2],Result[2],Result[0]]), np.append(CDF[MinInd:MaxInd+1],[CDF_Interpolate[2],0,0]), color=[0.8,0.8,0.8])
plt.plot([Result[2],Result[2]],[CDF_Interpolate[2],0],'k--')
plt.plot([Result[0],Result[0]],[CDF_Interpolate[0],0],'k--')
plt.plot([Result[1],Result[1]],[CDF_Interpolate[1],0],'k-')
plt.xlim(20,100)
plt.ylim(0,1)
plt.xlabel('Beach Width (m)')
plt.ylabel('Cumulative Probability Density')

plt.tight_layout()
plt.savefig('./BH_linear_chain_PDFs.pdf')
plt.savefig('./BH_linear_chain_PDFs.png')


ind = np.argsort(NewLike)
print "Most Likely RR1", RetreatRate1[ind[-3]], "RR2", RetreatRate2[ind[-3]], "BW", BeachWidth[ind[-3]]
print "Likelihood", NewLike[ind[-3]]
print "-Log Likelihood", -np.log(NewLike[ind[-3]])
