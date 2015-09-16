# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 09:42:30 2015

@author: mhurst
"""

#IMPORT MODULES
import numpy as np
import matplotlib.pyplot as plt

#setup file name for HG

Filename = "HG_change_chain2.dout"
print Filename

#Load the chain file
i, RetreatRate1, RetreatRate2, ChangeTime, BeachWidth, Elev_Init, NewLike, LastLike, NAccepted, NRejected  = np.loadtxt(Filename, unpack=True, skiprows=2, usecols=(0,1,2,3,4,5,6,7,8,9))

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

#Set number of bins
Numbins = 100
CDF_Interpolate = [0.025,0.5,0.975]

### PROCESS AND PLOT EVERYTHING ###

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

## function defining scientific notation
def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

#create dummy for blank axis labels
blanklabels = []

# Now plot the progression of parameters through time
# plot the number of accepted and reject parameters through time too
fig1 = plt.figure(1,figsize=(7.5,9.0))
plt.subplots_adjust(0.1,0.1,0.9,0.9)
ax = plt.subplot(511)
plt.plot(i,RetreatRate1,'k-')
ax.get_xaxis().set_ticks([])
plt.ylim(0.05,0.1)
plt.ylabel('Retreat Rate 1 (m/y)')
ax = plt.subplot(512)
plt.plot(i,RetreatRate2,'k-')
ax.get_xaxis().set_ticks([])
plt.ylabel('Retreat Rate 2 (m/y)')
ax = plt.subplot(513)
plt.plot(i,ChangeTime,'k-')
ax.get_xaxis().set_ticks([])
plt.ylabel('Change Time (years BP)')
ax = plt.subplot(514)
plt.plot(i,BeachWidth,'k-')
ax.get_xaxis().set_ticks([])
plt.ylabel('Beach Width (m)')
ax = plt.subplot(515)
plt.plot(i,np.abs(np.log(NewLike)),'k-')
plt.xlabel('Number of Iterations')
plt.ylabel('Negative Log-Likelihood')
plt.suptitle('All Results')
plt.tight_layout()
plt.savefig("./HG_change_chain_all.pdf")
plt.savefig("./HG_change_chain_all.png")


# Now plot the progression of parameters through time
# plot the number of accepted and reject parameters through time too
fig2 = plt.figure(2,figsize=(7.5,9.))
plt.subplots_adjust(0.1,0.1,0.9,0.9)

ax = plt.subplot(511)
plt.plot(Acceptedi,AcceptedRetreatRate1,'k-')
ax.get_xaxis().set_ticks([])
#plt.ylim(0.02,0.2)
plt.ylabel('Retreat Rate 1 (m/y)')

plt.axes([0.15,0.92,0.2,0.05])
plt.plot(Acceptedi[0:1000],AcceptedRetreatRate1[0:1000],'k-')
plt.xlim(0,1000)
plt.yticks([0.04,0.06,0.08,0.1])

ax = plt.subplot(512)
plt.plot(Acceptedi,AcceptedRetreatRate2,'k-')
ax.xaxis.set_ticklabels([])
#plt.ylim(0.,0.4)
plt.ylabel('Retreat Rate 2 (m/y)')

plt.axes([0.15,0.73,0.2,0.05])
plt.plot(Acceptedi[0:1000],AcceptedRetreatRate2[0:1000],'k-')
plt.xlim(0,1000)
plt.yticks([0.00,0.02,0.04,0.06,0.08,0.1])

ax = plt.subplot(513)
plt.plot(Acceptedi,AcceptedChangeTime,'k-')
ax.xaxis.set_ticklabels([])
plt.ylim(200,700)
plt.ylabel('Change Time (years BP)')

plt.axes([0.15,0.54,0.2,0.05])
plt.plot(Acceptedi[0:1000],AcceptedChangeTime[0:1000],'k-')
plt.xlim(0,1000)
plt.yticks([300,400,500,600])

ax = plt.subplot(514)
plt.plot(Acceptedi,AcceptedBeachWidth,'k-')
ax.xaxis.set_ticklabels([])
plt.ylabel('Beach Width (m)')

plt.axes([0.15,0.28,0.2,0.05])
plt.plot(Acceptedi[0:1000],AcceptedBeachWidth[0:1000],'k-')
plt.xlim(0,1000)
plt.yticks([0,20,40,60])

ax = plt.subplot(515)
plt.plot(i,np.abs(np.log(LastLike)),'k-')
plt.xlabel('Number of Iterations')
plt.ylabel('Negative Log-Likelihood')

plt.axes([0.15,0.15,0.2,0.05])
plt.plot(i[0:1000],np.abs(np.log(LastLike[0:1000])),'k-')
plt.xlim(0,1000)

plt.suptitle('Accepted Results')
plt.tight_layout()
plt.savefig("./HG_change_chain_accepted.pdf")
plt.savefig("./HG_change_chain_accepted.png")

########## DO PDF AND CDF plots ############

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

fig3 = plt.figure(3,figsize=(7.5,9))
plt.subplot(421)
plt.bar(Edges[:-1],Values,width = Edges[1]-Edges[0])
plt.xlim(0.04,0.08)
plt.xlabel('Retreat Rate 1 (m/y)')
plt.ylabel('Probability Density')

plt.subplot(422)
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
plt.xlim(0.04,0.08)
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

plt.subplot(423)
plt.bar(Edges[:-1],Values,width = Edges[1]-Edges[0])
plt.xlim(0.01,0.08)
plt.xlabel('Retreat Rate 2 (m/y)')
plt.ylabel('Probability Density')

plt.subplot(424)
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
plt.xlim(0.01,0.08)
plt.ylim(0,1)
plt.xlabel('Retreat Rate 2 (m/y)')
plt.ylabel('Cumulative Probability Density')

#Now for Change Time
Values, Edges = np.histogram(ChangeTime,30,normed=True,weights=NewLike)
Values /= np.sum(Values)
CDF = np.cumsum(Values)
LeftEdges = Edges[:-1]
RightEdges = Edges[1:]
Midpoints = (RightEdges + LeftEdges)/2
Result = np.interp(CDF_Interpolate,CDF,Midpoints)
print "Change Time:", Result
LikelyCT = Result[1]

plt.subplot(425)
plt.bar(Edges[:-1],Values,width = Edges[1]-Edges[0])
plt.xlim(0,1000)
plt.xlabel('Change Time (years BP)')
plt.ylabel('Probability Density')

plt.subplot(426)
plt.plot(Midpoints,CDF,'k-')
MinInd = np.argmin(np.abs(CDF-0.025))
MaxInd = np.argmin(np.abs(CDF-0.975))
MedInd = np.argmin(np.abs(CDF-0.5))
plt.plot([0,Result[0]],[CDF_Interpolate[0],CDF_Interpolate[0]],'r--')
plt.plot([0,Result[2]],[CDF_Interpolate[2],CDF_Interpolate[2]],'r--')
plt.plot([0,Result[1]],[CDF_Interpolate[1],CDF_Interpolate[1]],'r-')
plt.ylim(0,1)

plt.fill(np.append(Midpoints[MinInd:MaxInd+1],[Result[2],Result[2],Result[0]]), np.append(CDF[MinInd:MaxInd+1],[CDF_Interpolate[2],0,0]), color=[0.8,0.8,0.8])
plt.plot([Result[2],Result[2]],[CDF_Interpolate[2],0],'k--')
plt.plot([Result[0],Result[0]],[CDF_Interpolate[0],0],'k--')
plt.plot([Result[1],Result[1]],[CDF_Interpolate[1],0],'k-')
plt.xlim(0,1000)
plt.xlabel('Change Time (Years BP)')
plt.ylabel('Cumulative Probability Density')

#Finally for beach width
Values, Edges = np.histogram(BeachWidth,150,normed=True,weights=NewLike)
Values /= np.sum(Values)
CDF = np.cumsum(Values)
LeftEdges = Edges[:-1]
RightEdges = Edges[1:]
Midpoints = (RightEdges + LeftEdges)/2
Result = np.interp(CDF_Interpolate,CDF,Midpoints)
print "Beach Width:", Result
LikelyBW = Result[1]

plt.subplot(427)
plt.bar(Edges[:-1],Values,width = Edges[1]-Edges[0])
plt.xlim(35,60)
plt.xlabel('BeachWidth (m)')
plt.ylabel('Probability Density')

plt.subplot(428)
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
plt.xlim(35,60)
plt.ylim(0,1)
plt.xlabel('Beach Width (m)')
plt.ylabel('Cumulative Probability Density')

plt.tight_layout()
plt.savefig('./HG_change_chain_PDFs.pdf')
plt.savefig('./HG_change_chain_PDFs.png')


ind = np.argmax(NewLike)
print "Most Likely RR1", RetreatRate1[ind], "RR2", RetreatRate2[ind], "BW", BeachWidth[ind], "CT", ChangeTime[ind]
print "Likelihood", NewLike[ind]
print "-Log Likelihood", -np.log(NewLike[ind])

    
#####################################################
#       NOW PLOT change CHAIN FOR BEACHY HEAD       #
#####################################################

#setup file name for BH

Filename = "BH_change_chain2.dout"
print Filename

#Load the chain file
i, RetreatRate1, RetreatRate2, ChangeTime, BeachWidth, Elev_Init, NewLike, LastLike, NAccepted, NRejected  = np.loadtxt(Filename, unpack=True, skiprows=2, usecols=(0,1,2,3,4,5,6,7,8,9))

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

#Set number of bins
Numbins = 200
CDF_Interpolate = [0.025,0.5,0.975]

### PROCESS AND PLOT EVERYTHING ###

#create dummy for blank axis labels
blanklabels = []

# Now plot the progression of parameters through time
# plot the number of accepted and reject parameters through time too
fig4 = plt.figure(6,figsize=(7.5,9.0))
plt.subplots_adjust(0.1,0.1,0.9,0.9)
ax = plt.subplot(511)
plt.plot(i,RetreatRate1,'k-')
ax.xaxis.set_ticklabels([])
plt.ylabel('Retreat Rate 1 (m/y)')
ax = plt.subplot(512)
plt.plot(i,RetreatRate2,'k-')
ax.xaxis.set_ticklabels([])
plt.ylim(0.1,0.5)
plt.ylabel('Retreat Rate 2 (m/y)')
ax = plt.subplot(513)
plt.plot(i,ChangeTime,'k-')
ax.xaxis.set_ticklabels([])
plt.ylabel('Change Time (years BP)')
ax = plt.subplot(514)
plt.plot(i,BeachWidth,'k-')
ax.xaxis.set_ticklabels([])
plt.ylabel('Beach Width (m)')
plt.subplot(515)
plt.plot(i,np.abs(np.log(NewLike)),'k-')
plt.xlabel('Number of Iterations')
plt.ylabel('Negative Log-Likelihood')
plt.suptitle('All Results')
plt.tight_layout()
plt.savefig("./BH_change_chain_all.pdf")
plt.savefig("./BH_change_chain_all.png")

# Now plot the progression of parameters through time
# plot the number of accepted and reject parameters through time too
fig5 = plt.figure(7,figsize=(7.5,9.))
plt.subplots_adjust(0.1,0.1,0.9,0.9)

ax = plt.subplot(511)
plt.plot(Acceptedi,AcceptedRetreatRate1,'k-')
ax.get_xaxis().set_ticks([])
#plt.ylim(0.04,0.12)
plt.ylabel('Retreat Rate 1 (m/y)')

plt.axes([0.15,0.92,0.2,0.05])
plt.plot(Acceptedi[0:1000],AcceptedRetreatRate1[0:1000],'k-')
plt.xlim(0,1000)
plt.yticks([0.04,0.06,0.08,0.1])

ax = plt.subplot(512)
plt.plot(Acceptedi,AcceptedRetreatRate2,'k-')
ax.xaxis.set_ticklabels([])
plt.ylim(0.10,0.50)
plt.ylabel('Retreat Rate 2 (m/y)')

plt.axes([0.15,0.73,0.2,0.05])
plt.plot(Acceptedi[0:1000],AcceptedRetreatRate2[0:1000],'k-')
plt.xlim(0,1000)
#plt.yticks([0.04,0.06,0.08,0.1])

ax = plt.subplot(513)
plt.plot(Acceptedi,AcceptedChangeTime,'k-')
ax.xaxis.set_ticklabels([])
plt.ylim(200,1200)
plt.ylabel('Change Time (years BP)')

plt.axes([0.15,0.54,0.2,0.05])
plt.plot(Acceptedi[0:1000],AcceptedChangeTime[0:1000],'k-')
plt.xlim(0,1000)
plt.yticks([200,400,600,800])

ax = plt.subplot(514)
plt.plot(Acceptedi,AcceptedBeachWidth,'k-')
ax.xaxis.set_ticklabels([])
plt.ylim(0,100)
plt.ylabel('Beach Width (m)')

plt.axes([0.15,0.35,0.2,0.05])
plt.plot(Acceptedi[0:1000],AcceptedBeachWidth[0:1000],'k-')
plt.xlim(0,1000)
plt.yticks([20,40,60,80,100])

ax = plt.subplot(515)
plt.plot(i,np.abs(np.log(LastLike)),'k-')
plt.ylim(80,120)
plt.xlabel('Number of Iterations')
plt.ylabel('Negative Log-Likelihood')

plt.axes([0.15,0.16,0.2,0.05])
plt.plot(i[0:1000],np.abs(np.log(LastLike[0:1000])),'k-')
plt.xlim(0,1000)
plt.yticks([80,100,120,140])

plt.suptitle('Accepted Results')
plt.tight_layout()
plt.savefig("./BH_change_chain_accepted.pdf")
plt.savefig("./BH_change_chain_accepted.png")


########### DO PDF AND CDF plots ############

#First for Retreat Rate 1
Values, Edges = np.histogram(RetreatRate1,200,normed=True,weights=NewLike)
Values /= np.sum(Values)
CDF = np.cumsum(Values)
LeftEdges = Edges[:-1]
RightEdges = Edges[1:]
Midpoints = (RightEdges + LeftEdges)/2
Result = np.interp(CDF_Interpolate,CDF,Midpoints)
print "Retreat Rate 1:", Result
LikelyRR1 = Result[1]

fig3 = plt.figure(8,figsize=(7.5,9))
plt.subplot(421)
plt.bar(Edges[:-1],Values,width = Edges[1]-Edges[0])
plt.xlim(0.02,0.035)
plt.xlabel('Retreat Rate 1 (m/y)')
plt.ylabel('Probability Density')

plt.subplot(422)
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
plt.xlim(0.02,0.035)
plt.ylim(0,1)
plt.xlabel('Retreat Rate 1 (m/y)')
plt.ylabel('Cumulative Probability Density')

#Retreat Rate 2
Values, Edges = np.histogram(RetreatRate2,50,normed=True,weights=NewLike)
Values /= np.sum(Values)
CDF = np.cumsum(Values)
LeftEdges = Edges[:-1]
RightEdges = Edges[1:]
Midpoints = (RightEdges + LeftEdges)/2
Result = np.interp(CDF_Interpolate,CDF,Midpoints)
print "Retreat Rate 2:", Result
LikelyRR2 = Result[1]

plt.subplot(423)
plt.bar(Edges[:-1],Values,width = Edges[1]-Edges[0])
plt.xlim(0.15,0.45)
plt.xlabel('Retreat Rate 2 (m/y)')
plt.ylabel('Probability Density')

plt.subplot(424)
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
plt.xlim(0.15,0.45)
plt.ylim(0,1)
plt.xlabel('Retreat Rate 2 (m/y)')
plt.ylabel('Cumulative Probability Density')

#Now for Change Time
Values, Edges = np.histogram(ChangeTime,50,normed=True,weights=NewLike)
Values /= np.sum(Values)
CDF = np.cumsum(Values)
LeftEdges = Edges[:-1]
RightEdges = Edges[1:]
Midpoints = (RightEdges + LeftEdges)/2
Result = np.interp(CDF_Interpolate,CDF,Midpoints)
print "Change Time:", Result
LikelyCT = Result[1]

plt.subplot(425)
plt.bar(Edges[:-1],Values,width = Edges[1]-Edges[0])
plt.xlim(0,750)
plt.xlabel('Change Time (years BP)')
plt.ylabel('Probability Density')

plt.subplot(426)
plt.plot(Midpoints,CDF,'k-')
MinInd = np.argmin(np.abs(CDF-0.025))
MaxInd = np.argmin(np.abs(CDF-0.975))
MedInd = np.argmin(np.abs(CDF-0.5))
plt.plot([0,Result[0]],[CDF_Interpolate[0],CDF_Interpolate[0]],'r--')
plt.plot([0,Result[2]],[CDF_Interpolate[2],CDF_Interpolate[2]],'r--')
plt.plot([0,Result[1]],[CDF_Interpolate[1],CDF_Interpolate[1]],'r-')
plt.ylim(0,1)

plt.fill(np.append(Midpoints[MinInd:MaxInd+1],[Result[2],Result[2],Result[0]]), np.append(CDF[MinInd:MaxInd+1],[CDF_Interpolate[2],0,0]), color=[0.8,0.8,0.8])
plt.plot([Result[2],Result[2]],[CDF_Interpolate[2],0],'k--')
plt.plot([Result[0],Result[0]],[CDF_Interpolate[0],0],'k--')
plt.plot([Result[1],Result[1]],[CDF_Interpolate[1],0],'k-')
plt.xlim(0,750)
plt.xlabel('Change Time (Years BP)')
plt.ylabel('Cumulative Probability Density')

#Finally for beach width
Values, Edges = np.histogram(BeachWidth,50,normed=True,weights=NewLike)
Values /= np.sum(Values)
CDF = np.cumsum(Values)
LeftEdges = Edges[:-1]
RightEdges = Edges[1:]
Midpoints = (RightEdges + LeftEdges)/2
Result = np.interp(CDF_Interpolate,CDF,Midpoints)
LikelyBW = Result[1]
print "Beach Width:", Result

plt.subplot(427)
plt.bar(Edges[:-1],Values,width = Edges[1]-Edges[0])
plt.xlim(0,40)
plt.xlabel('Beach Width (m)')
plt.ylabel('Probability Density')

plt.subplot(428)
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
plt.xlim(0,40)
plt.xlabel('Beach Width (m)')
plt.ylabel('Cumulative Probability Density')
plt.ylim(0,1)
plt.tight_layout()
plt.savefig('./BH_change_chain_PDFs.pdf')
plt.savefig('./BH_change_chain_PDFs.png')


ind = np.argmax(NewLike)
print "Most Likely RR1", RetreatRate1[ind], "RR2", RetreatRate2[ind], "BW", BeachWidth[ind]
print "Likelihood", NewLike[ind]
print "-Log Likelihood", -np.log(NewLike[ind])

###
# PLOT RR2 vs CT
###

plt.figure(99,figsize=(3.5,3.5))
plt.plot(AcceptedRetreatRate2[::2],AcceptedChangeTime[::2],'k.')
plt.xlabel(r'Retreat Rate 2 (m yr\textsuperscript{-1})')
plt.ylabel('Change Time (years BP)')
plt.tight_layout()
plt.savefig("BH_change_RR2_vs_CT.pdf")
plt.savefig("BH_change_RR2_vs_CT.png")
#######################################
###
### Dataset has two distinct peaks, split into two based on 0.06m/y threshold
### and analyse PDFs and CDFs separately
###
#######################################
#
#RR1_Peak1 = RetreatRate1[RetreatRate1 < 0.013]
#RR2_Peak1 = RetreatRate2[RetreatRate1 < 0.013]
#CT_Peak1 = ChangeTime[RetreatRate1 < 0.013]
#BW_Peak1 = BeachWidth[RetreatRate1 < 0.013]
#NL_Peak1 = NewLike[RetreatRate1 < 0.013]
#
#Values, Edges = np.histogram(RR1_Peak1,40,normed=True,weights=NL_Peak1)
#Values /= np.sum(Values)
#CDF = np.cumsum(Values)
#LeftEdges = Edges[:-1]
#RightEdges = Edges[1:]
#Midpoints = (RightEdges + LeftEdges)/2
#Result = np.interp(CDF_Interpolate,CDF,Midpoints)
#print "Retreat Rate BH Peak 1:", Result
#
#ind = np.argmax(NL_Peak1)
#print "Most Likely RR", RR1_Peak1[ind], "BW", BW_Peak1[ind]
#print "Likelihood", NL_Peak1[ind]
#print "-Log Likelihood", -np.log(NL_Peak1[ind])
#
#plt.figure(4,figsize=(7.5,9))
#plt.subplot(421)
#plt.bar(Edges[:-1],Values,width = Edges[1]-Edges[0])
#plt.xlim(0.010,0.014)
#plt.xticks([0.010,0.011,0.012,0.013,0.014])
#plt.xlabel('Retreat Rate (m/y)')
#plt.ylabel('Probability Density')
#
#plt.subplot(422)
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
#plt.xlim(0.01,0.014)
#plt.ylim(0,1)
#plt.xlabel('Retreat Rate (m/y)')
#plt.ylabel('Cumulative Probability Density')
#
#Values, Edges = np.histogram(RR2_Peak1,100,normed=True,weights=NL_Peak1)
#Values /= np.sum(Values)
#CDF = np.cumsum(Values)
#LeftEdges = Edges[:-1]
#RightEdges = Edges[1:]
#Midpoints = (RightEdges + LeftEdges)/2
#Result = np.interp(CDF_Interpolate,CDF,Midpoints)
#print "Retreat Rate 2:", Result
#
#plt.subplot(423)
#plt.bar(Edges[:-1],Values,width = Edges[1]-Edges[0])
#plt.xlim(0.2,0.3)
#plt.xlabel('Retreat Rate 2 (m/yr)')
#plt.ylabel('Probability Density')
#
#plt.subplot(424)
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
#plt.xlim(0.2,0.3)
#plt.ylim(0,1)
#plt.xlabel('Retreat Rate 2 (m/y)')
#plt.ylabel('Cumulative Probability Density')
#
#Values, Edges = np.histogram(CT_Peak1,100,normed=True,weights=NL_Peak1)
#Values /= np.sum(Values)
#CDF = np.cumsum(Values)
#LeftEdges = Edges[:-1]
#RightEdges = Edges[1:]
#Midpoints = (RightEdges + LeftEdges)/2
#Result = np.interp(CDF_Interpolate,CDF,Midpoints)
#print "Change Time:", Result
#
#plt.subplot(425)
#plt.bar(Edges[:-1],Values,width = Edges[1]-Edges[0])
#plt.xlim(350,700)
#plt.xlabel('Change Time (y BP)')
#plt.ylabel('Probability Density')
#
#plt.subplot(426)
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
#plt.xlim(350,700)
#plt.ylim(0,1)
#plt.xlabel('Change Time (y BP)')
#plt.ylabel('Cumulative Probability Density')
#
#Values, Edges = np.histogram(BW_Peak1,200,normed=True,weights=NL_Peak1)
#Values /= np.sum(Values)
#CDF = np.cumsum(Values)
#LeftEdges = Edges[:-1]
#RightEdges = Edges[1:]
#Midpoints = (RightEdges + LeftEdges)/2
#Result = np.interp(CDF_Interpolate,CDF,Midpoints)
#print "Beach Width:", Result
#
#plt.subplot(427)
#plt.bar(Edges[:-1],Values,width = Edges[1]-Edges[0])
#plt.xlim(0,10)
#plt.xlabel('BeachWidth (m)')
#plt.ylabel('Probability Density')
#
#plt.subplot(428)
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
#plt.xlim(0,10)
#plt.ylim(0,1)
#plt.xlabel('Beach Width (m)')
#plt.ylabel('Cumulative Probability Density')
#
#plt.tight_layout()
#plt.savefig('./BH_change_chain_PDFs_Peak1.pdf')
#plt.savefig('./BH_change_chain_PDFs_Peak1.png')



########## NOW FOR PEAK 2 #############
#
#
#indicies = np.where(RetreatRate1>0.013)
#RR1_Peak2 = RetreatRate1[indicies]
#RR2_Peak2 = RetreatRate2[indicies]
#CT_Peak2 = ChangeTime[indicies]
#BW_Peak2 = BeachWidth[indicies]
#NL_Peak2 = NewLike[indicies]
#
#Values, Edges = np.histogram(RR1_Peak2,200,normed=True,weights=NL_Peak2)
#Values /= np.sum(Values)
#CDF = np.cumsum(Values)
#LeftEdges = Edges[:-1]
#RightEdges = Edges[1:]
#Midpoints = (RightEdges + LeftEdges)/2
#Result = np.interp(CDF_Interpolate,CDF,Midpoints)
#print "Retreat Rate BH Peak 2:", Result
#
#ind = np.argmax(NL_Peak2)
#print "Most Likely RR", RR1_Peak2[ind], "BW", BW_Peak2[ind]
#print "Likelihood", NL_Peak2[ind]
#print "-Log Likelihood", -np.log(NL_Peak2[ind])
#
#plt.figure(5,figsize=(7.5,9))
#plt.subplot(421)
#plt.bar(Edges[:-1],Values,width = Edges[1]-Edges[0])
#plt.xlim(0.012,0.035)
#plt.xlabel('Retreat Rate 1 (m/y)')
#plt.ylabel('Probability Density')
#
#plt.subplot(422)
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
#plt.xlim(0.012,0.035)
#plt.ylim(0,1)
#plt.xlabel('Retreat Rate 1 (m/y)')
#plt.ylabel('Cumulative Probability Density')
#
#Values, Edges = np.histogram(RR2_Peak2,200,normed=True,weights=NL_Peak2)
#Values /= np.sum(Values)
#CDF = np.cumsum(Values)
#LeftEdges = Edges[:-1]
#RightEdges = Edges[1:]
#Midpoints = (RightEdges + LeftEdges)/2
#Result = np.interp(CDF_Interpolate,CDF,Midpoints)
#print "Retreat Rate 2:", Result
#
#plt.subplot(423)
#plt.bar(Edges[:-1],Values,width = Edges[1]-Edges[0])
#plt.xlim(0.2,0.4)
#
#plt.xlabel('Retreat Rate 2 (m/y)')
#plt.ylabel('Probability Density')
#
#plt.subplot(424)
#plt.plot(Edges[:-1],CDF,'k-')
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
#plt.ylim(0,1)
#plt.xlim(0.2,0.4)
#plt.xlabel('Retreat Rate 2 (m/y)')
#plt.ylabel('Cumulative Probability Density')
#
#
#Values, Edges = np.histogram(CT_Peak2,100,normed=True,weights=NL_Peak2)
#Values /= np.sum(Values)
#CDF = np.cumsum(Values)
#LeftEdges = Edges[:-1]
#RightEdges = Edges[1:]
#Midpoints = (RightEdges + LeftEdges)/2
#Result = np.interp(CDF_Interpolate,CDF,Midpoints)
#print "Change Time:", Result
#
#plt.subplot(425)
#plt.bar(Edges[:-1],Values,width = Edges[1]-Edges[0])
#plt.xlim(200.,800)
#
#plt.xlabel('Change Time (y BP)')
#plt.ylabel('Probability Density')
#
#plt.subplot(426)
#plt.plot(Edges[:-1],CDF,'k-')
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
#plt.ylim(0,1)
#plt.xlim(200,800)
#plt.xlabel('Change Time (y BP)')
#plt.ylabel('Cumulative Probability Density')
#
#Values, Edges = np.histogram(BW_Peak2,Numbins,normed=True,weights=NL_Peak2)
#Values /= np.sum(Values)
#CDF = np.cumsum(Values)
#LeftEdges = Edges[:-1]
#RightEdges = Edges[1:]
#Midpoints = (RightEdges + LeftEdges)/2
#Result = np.interp(CDF_Interpolate,CDF,Midpoints)
#print "Beach Width:", Result
#
#plt.subplot(427)
#plt.bar(Edges[:-1],Values,width = Edges[1]-Edges[0])
#plt.xlim(0.,30)
#
#plt.xlabel('BeachWidth (m)')
#plt.ylabel('Probability Density')
#
#plt.subplot(428)
#plt.plot(Edges[:-1],CDF,'k-')
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
#plt.ylim(0,1)
#plt.xlim(0,30)
#plt.xlabel('Beach Width (m)')
#plt.ylabel('Cumulative Probability Density')
#
#plt.tight_layout()
#plt.savefig('./BH_change_chain_PDFs_Peak2.pdf')
#plt.savefig('./BH_change_chain_PDFs_Peak2.png')


