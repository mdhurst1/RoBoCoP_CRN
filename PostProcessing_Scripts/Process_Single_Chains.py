# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 09:42:30 2015

@author: mhurst
"""

#IMPORT MODULES
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


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


###########################
### START WITH HOPE GAP ###
###########################

#setup file name for HG    
Filename = "HG_single_chain2.dout"
print "Processing MCMC Chain Output:", Filename

#Load the chain file
i, RetreatRate1, BeachWidth, NewLike, LastLike, NAccepted, NRejected  = np.loadtxt(Filename, unpack=True, skiprows=2, usecols=(0,1,4,6,7,8,9))

# Get accepted parameters for plotting convergence
N = 0
AcceptedMask = np.zeros(len(i),dtype=np.int)
for j in range(0,len(i)):
    if NAccepted[j] > N:
        N += 1
        AcceptedMask[j] = 1

AcceptedRetreatRate1 = RetreatRate1[AcceptedMask == 1]
AcceptedBeachWidth = BeachWidth[AcceptedMask == 1]
Acceptedi = i[AcceptedMask == 1]

#Set number of bins
Numbins = 100
CDF_Interpolate = [0.025,0.5,0.975]
print "Percentiles", CDF_Interpolate

# Now plot the progression of parameters through time
# plot the number of accepted and reject parameters through time too
fig1 = plt.figure(1,figsize=(7.5,5.125))
plt.subplots_adjust(0.1,0.1,0.9,0.9)
plt.subplot(311)
plt.plot(i,RetreatRate1,'k-')
plt.xlabel('Number of Iterations')
plt.ylabel('Retreat Rate 1 (m/y)')
plt.subplot(312)
plt.plot(i,BeachWidth,'k-')
plt.xlabel('Number of Iterations')
plt.ylabel('Beach Width (m)')
plt.subplot(313)
plt.plot(i,-np.log(NewLike),'k-')
plt.xlabel('Number of Iterations')
plt.ylabel('Negative Log-Likelihood')
plt.tight_layout()
plt.savefig("./HG_single_chain_all.pdf")
plt.savefig("./HG_single_chain_all.png")

# Now plot the progression of parameters through time
# plot the number of accepted and reject parameters through time too
fig2 = plt.figure(2,figsize=(7.5,5.125))
plt.subplots_adjust(0.1,0.1,0.9,0.9)

plt.subplot(311)
plt.plot(Acceptedi,AcceptedRetreatRate1,'k-')
plt.xlabel('Number of Iterations')
plt.ylabel('Retreat Rate 1 (m/y)')

plt.subplot(312)
plt.plot(Acceptedi,AcceptedBeachWidth,'k-')
plt.xlabel('Number of Iterations')
plt.ylabel('Beach Width (m)')

plt.subplot(313)
plt.plot(i,-np.log(LastLike),'k-')
#plt.ylim(10**(-30),10**(-14))
plt.xlabel('Number of Iterations')
plt.ylabel('Likelihood')

plt.axes([0.74,0.85,0.2,0.1])
plt.plot(Acceptedi[0:30],AcceptedRetreatRate1[0:30],'k-')
plt.xlim(0,100)
plt.yticks([0.07,0.08,0.09,0.1])

plt.axes([0.74,0.46,0.2,0.1])
plt.plot(Acceptedi[0:30],AcceptedBeachWidth[0:30],'k-')
plt.xlim(0,60)
plt.yticks([0,20,40,60])

plt.axes([0.74,0.2,0.2,0.1])
plt.plot(i[0:100],-np.log(LastLike[0:100]),'k-')
plt.yticks([40,60,80,100])
plt.ylim(40,100)
#plt.yticks([10**(-200),10**(-100),10**(0)])

plt.tight_layout()
plt.savefig("./HG_single_chain_accepted.pdf")
plt.savefig("./HG_single_chain_accepted.png")

### GOOD TO HERE ###

########### DO PDF AND CDF plots ############
Values, Edges = np.histogram(RetreatRate1,Numbins,normed=True,weights=NewLike)
Values /= np.sum(Values)
CDF = np.cumsum(Values)
LeftEdges = Edges[:-1]
RightEdges = Edges[1:]
Midpoints = (RightEdges + LeftEdges)/2
Result = np.interp(CDF_Interpolate,CDF,Midpoints)
print "Retreat Rate:", Result
LikelyRR1 = Result[1]

fig3 = plt.figure(3,figsize=(7.5,5.5))
plt.subplot(221)
plt.bar(Edges[:-1],Values,width = Edges[1]-Edges[0])
plt.xlim(0.04,0.08)
plt.xlabel('Retreat Rate (m/y)')
plt.ylabel('Probability Density')

plt.subplot(223)
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
plt.xlabel('Retreat Rate (m/y)')
plt.ylabel('Cumulative Probability Density')

Values, Edges = np.histogram(BeachWidth,Numbins,normed=True,weights=NewLike)
Values /= np.sum(Values)
CDF = np.cumsum(Values)
LeftEdges = Edges[:-1]
RightEdges = Edges[1:]
Midpoints = (RightEdges + LeftEdges)/2
Result = np.interp(CDF_Interpolate,CDF,Midpoints)
print "Beach Width:", Result
LikelyBW = Result[1]

plt.subplot(222)
plt.bar(Edges[:-1],Values,width = Edges[1]-Edges[0])
#plt.xlim(20,60)
plt.xlabel('BeachWidth (m)')
plt.ylabel('Probability Density')

plt.subplot(224)
plt.plot(Edges[:-1],CDF,'k-')
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
#plt.xlim(20,60)
plt.xlabel('Beach Width (m)')
plt.ylabel('Cumulative Probability Density')
plt.ylim(0,1)

plt.tight_layout()
plt.savefig('./HG_single_chain_PDFs.pdf')
plt.savefig('./HG_single_chain_PDFs.png')

#print most likely
ind = np.argmax(NewLike)
print "Most Likely RR1", RetreatRate1[ind],"BW", BeachWidth[ind]
print "Likelihood", NewLike[ind]
print "-Log Likelihood", -np.log(NewLike[ind])


######################################################
##       NOW PLOT SINGLE CHAIN FOR BEACHY HEAD       #
######################################################
#
Filename = "BH_single_chain2.dout"
print "Processing MCMC chain:", Filename

#Load the chain file
i, RetreatRate1, BeachWidth, NewLike, LastLike, NAccepted, NRejected  = np.loadtxt(Filename, unpack=True, skiprows=2, usecols=(0,1,4,6,7,8,9))
Numbins = 100
CDF_Interpolate = [0.025,0.5,0.975]

# Get accepted parameters for plotting convergence
N = 0
AcceptedMask = np.zeros(len(i),dtype=np.int)
for j in range(0,len(i)):
    if NAccepted[j] > N:
        N += 1
        AcceptedMask[j] = 1

AcceptedRetreatRate1 = RetreatRate1[AcceptedMask == 1]
AcceptedBeachWidth = BeachWidth[AcceptedMask == 1]
Acceptedi = i[AcceptedMask == 1]

# Now plot the progression of parameters through time
# plot the number of accepted and reject parameters through time too
fig5 = plt.figure(8,figsize=(7.5,5.125))
plt.subplots_adjust(0.1,0.1,0.9,0.9)
plt.subplot(311)
plt.plot(i,RetreatRate1,'k-')
plt.xlabel('Number of Iterations')
plt.ylabel('Retreat Rate 1 (m/y)')
plt.subplot(312)
plt.plot(i,BeachWidth,'k-')
plt.xlabel('Number of Iterations')
plt.ylabel('Beach Width (m)')
plt.subplot(313)
plt.plot(i,-np.log(NewLike),'k-')
plt.xlabel('Number of Iterations')
plt.ylabel('Negative Log-Likelihood')
plt.suptitle('All Results')
plt.tight_layout()
plt.savefig("./BH_single_chain_all.pdf")
plt.savefig("./BH_single_chain_all.png")

# Now plot the progression of parameters through time
# plot the number of accepted and reject parameters through time too
fig6 = plt.figure(9,figsize=(7.5,5.125))
plt.subplots_adjust(0.1,0.1,0.9,0.9)

plt.subplot(311)
plt.plot(Acceptedi,AcceptedRetreatRate1,'k-')
plt.xlabel('Number of Iterations')
plt.ylabel('Retreat Rate 1 (m/y)')
plt.ylim(0.04,0.1)

plt.subplot(312)
plt.plot(Acceptedi,AcceptedBeachWidth,'k-')
plt.xlabel('Number of Iterations')
plt.ylabel('Beach Width (m)')
plt.ylim(0.,50)

plt.subplot(313)
plt.plot(i,-np.log(LastLike),'k-')
plt.ylim(120,150)
plt.xlabel('Number of Iterations')
plt.ylabel('Negative Log-Likelihood')

plt.axes([0.18,0.85,0.2,0.1])
plt.plot(Acceptedi[0:1000],AcceptedRetreatRate1[0:1000],'k-')
plt.xlim(0,1000)
plt.yticks([0.04,0.06,0.08,0.10,0.12])

plt.axes([0.18,0.45,0.2,0.1])
plt.plot(Acceptedi[0:1000],AcceptedBeachWidth[0:1000],'k-')
plt.xlim(0,1000)
plt.yticks([0.,20.,40.,60.])

plt.axes([0.18,0.2,0.2,0.1])
plt.plot(i[0:1000],-np.log(LastLike[0:1000]),'k-')
plt.xlim(0,1000)
plt.ylim(120,150)
plt.yticks([120,130,140,150])

plt.suptitle('Accepted Results')
plt.tight_layout()
plt.savefig("./BH_single_chain_accepted.pdf")
plt.savefig("./BH_single_chain_accepted.png")

###############################
##        CDF FOR BH         ##
###############################

# FIRST FOR RETREAT RATE

RRValues, RREdges = np.histogram(RetreatRate1,Numbins,normed=True,weights=NewLike)
RRValues /= np.sum(RRValues)
RRCDF = np.cumsum(RRValues)
RRLeftEdges = RREdges[:-1]
RRRightEdges = RREdges[1:]
RRMidpoints = (RRRightEdges + RRLeftEdges)/2
RRResult = np.interp(CDF_Interpolate,RRCDF,RRMidpoints)
print "Retreat Rate:", RRResult

fig7 = plt.figure(10,figsize=(7.5,5.5))
plt.subplot(221)
plt.bar(RREdges[:-1],RRValues,width = RREdges[1]-RREdges[0])
plt.xlim(0.02,0.08)
plt.xlabel('Retreat Rate (m/y)')
plt.ylabel('Probability Density')

plt.subplot(223)
plt.plot(RRMidpoints,RRCDF,'k-')
MinInd = np.argmin(np.abs(RRCDF-0.025))
MaxInd = np.argmin(np.abs(RRCDF-0.975))
MedInd = np.argmin(np.abs(RRCDF-0.5))
plt.plot([0,RRResult[0]],[CDF_Interpolate[0],CDF_Interpolate[0]],'r--')
plt.plot([0,RRResult[2]],[CDF_Interpolate[2],CDF_Interpolate[2]],'r--')
plt.plot([0,RRResult[1]],[CDF_Interpolate[1],CDF_Interpolate[1]],'r-')

plt.fill(np.append(RRMidpoints[MinInd:MaxInd+1],[RRResult[2],RRResult[2],RRResult[0]]), np.append(RRCDF[MinInd:MaxInd+1],[CDF_Interpolate[2],0,0]), color=[0.8,0.8,0.8])
plt.plot([RRResult[2],RRResult[2]],[CDF_Interpolate[2],0],'k--')
plt.plot([RRResult[0],RRResult[0]],[CDF_Interpolate[0],0],'k--')
plt.plot([RRResult[1],RRResult[1]],[CDF_Interpolate[1],0],'k-')
plt.xlim(0.02,0.08)
plt.ylim(0,1)
plt.xlabel('Retreat Rate (m/y)')
plt.ylabel('Cumulative Probability Density')


# NOW FOR BEACH WIDTH

WidthValues, WidthEdges = np.histogram(BeachWidth,Numbins,normed=True,weights=NewLike)
WidthValues /= np.sum(WidthValues)
WidthCDF = np.cumsum(WidthValues)
WidthLeftEdges = WidthEdges[:-1]
WidthRightEdges = WidthEdges[1:]
WidthMidpoints = (WidthRightEdges + WidthLeftEdges)/2
WidthResult = np.interp(CDF_Interpolate,WidthCDF,WidthMidpoints)
print "Beach Width:", WidthResult
ind = np.argmax(NewLike)
print "most Likely RR", RetreatRate1[ind], "BW", BeachWidth[ind]
print "Likelihood", NewLike[ind]
print "-Log Likelihood", -np.log(NewLike[ind])

plt.subplot(222)
plt.bar(WidthEdges[:-1],WidthValues,width = WidthEdges[1]-WidthEdges[0])
plt.xlim(20,80)
plt.xlabel('BeachWidth (m)')
plt.ylabel('Probability Density')

plt.subplot(224)
plt.plot(WidthEdges[:-1],WidthCDF,'k-')
MinInd = np.argmin(np.abs(WidthCDF-0.025))
MaxInd = np.argmin(np.abs(WidthCDF-0.975))
MedInd = np.argmin(np.abs(WidthCDF-0.5))
plt.plot([0,WidthResult[0]],[CDF_Interpolate[0],CDF_Interpolate[0]],'r--')
plt.plot([0,WidthResult[2]],[CDF_Interpolate[2],CDF_Interpolate[2]],'r--')
plt.plot([0,WidthResult[1]],[CDF_Interpolate[1],CDF_Interpolate[1]],'r-')

plt.fill(np.append(WidthMidpoints[MinInd:MaxInd+1],[WidthResult[2],WidthResult[2],WidthResult[0]]), np.append(WidthCDF[MinInd:MaxInd+1],[CDF_Interpolate[2],0,0]), color=[0.8,0.8,0.8])
plt.plot([WidthResult[2],WidthResult[2]],[CDF_Interpolate[2],0],'k--')
plt.plot([WidthResult[0],WidthResult[0]],[CDF_Interpolate[0],0],'k--')
plt.plot([WidthResult[1],WidthResult[1]],[CDF_Interpolate[1],0],'k-')
plt.xlim(20,80)
plt.ylim(0,1)
plt.xlabel('Beach Width (m)')
plt.ylabel('Cumulative Probability Density')

plt.tight_layout()
plt.savefig('./BH_single_chain_PDFs.pdf')
plt.savefig('./BH_single_chain_PDFs.png')

#######################################
## plot 2d histogram contour plot    ##
#######################################

Values, Xedges, Yedges = np.histogram2d(RetreatRate1,BeachWidth,bins=400,normed=True,weights=NewLike)
Values /= np.sum(Values)

XLeftEdges = Xedges[:-1]
XRightEdges = Xedges[1:]
XMidpoints = (XRightEdges+XLeftEdges)/2

YLeftEdges = Yedges[:-1]
YRightEdges = Yedges[1:]
YMidpoints = (YRightEdges+YLeftEdges)/2 
X,Y = np.meshgrid(XMidpoints,YMidpoints)

def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)
    
plt.figure(11,figsize=(3.75,2.75))
extent = [Xedges[0], Xedges[-1], Yedges[0], Yedges[-1]]
plt.imshow(Values.T,cmap='Reds',extent=extent,origin='lower',aspect='auto')
plt.xlim(0.02,0.1)
plt.ylim(10,60)
plt.xlabel("Retreat Rate (m/y)")
plt.ylabel("Beach Width (m)")
plt.colorbar(format=ticker.FuncFormatter(fmt),ticks=[0,0.002,0.004,0.006],label='Probability Density')


plt.tight_layout()
plt.savefig('./BH_2DHist.pdf')
plt.savefig('./BH_2DHist.png')
#plt.xlim(30,80)
#plt.ylim(0)

