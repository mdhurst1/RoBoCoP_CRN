# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 17:22:23 2016

Tidal constituents

@author: mhurst
"""


#import modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

# Customise figure style #
rc('font',size=8)
rc('ytick.major',pad=5)
rc('xtick.major',pad=5)
rc('text', usetex=True)

#The equillibrium tide
M1_Speed = 14.492
M1_Amplitude = 0.1
M2_Speed = 28.984
M2_Amplitude = 2.5 
O1_Speed = 13.943
O1_Amplitude = 0.5
K1_Speed = 15.041
K1_Amplitude = 0.4
N2_Amplitude = 0.5
N2_Speed = 28.440
K2_Speed = 30.082
K2_Amplitude = 0.1

#convert dates and times to datetimes
TideTimes = np.arange('2016-01-01T00:00','2016-02-01T01:00',dtype='datetime64[h]')
TimesHours = np.arange(0,len(TideTimes))

##Tides
WaterLevels = np.zeros(len(TideTimes))

#Setup figure
plt.figure(1,figsize=(6.6,6.6))

#simple diurnal
K1_Amplitude = 2.5
WaterLevels = K1_Amplitude*np.cos(np.radians(K1_Speed*TimesHours))
ax = plt.subplot(611)
plt.plot(TideTimes,WaterLevels,'k-')
ax.set_xticklabels([])
plt.yticks([-2,-1,0,1,2])
plt.ylabel("Sea Level (m)")

##Simple semi-diurnal
WaterLevels = M2_Amplitude*np.cos(np.radians(M2_Speed*TimesHours))
ax = plt.subplot(612)
plt.plot(TideTimes,WaterLevels,'k-')
ax.set_xticklabels([])
plt.yticks([-2,-1,0,1,2])
plt.ylabel("Sea Level (m)")

##Simple mixed
M2_Amplitude = 1.8
K1_Amplitude = 0.7
WaterLevels = K1_Amplitude*np.cos(np.radians(K1_Speed*TimesHours))
WaterLevels += M2_Amplitude*np.cos(np.radians(M2_Speed*TimesHours))
ax = plt.subplot(613)
plt.plot(TideTimes,WaterLevels,'k-')
ax.set_xticklabels([])
plt.yticks([-2,-1,0,1,2])
plt.ylabel("Sea Level (m)")

#Lunar diurnal
K1_Amplitude = 2.
WaterLevels = K1_Amplitude*np.cos(np.radians(K1_Speed*TimesHours))
WaterLevels += O1_Amplitude*np.cos(np.radians(O1_Speed*TimesHours))
ax = plt.subplot(614)
plt.plot(TideTimes,WaterLevels,'k-')
ax.set_xticklabels([])
plt.yticks([-2,-1,0,1,2])
plt.ylabel("Sea Level (m)")

##lunar semidiurnal
M2_Amplitude = 2.
WaterLevels = M2_Amplitude*np.cos(np.radians(M2_Speed*TimesHours))
WaterLevels += N2_Amplitude*np.cos(np.radians(N2_Speed*TimesHours))
ax = plt.subplot(615)
plt.plot(TideTimes,WaterLevels,'k-')
ax.set_xticklabels([])
plt.yticks([-2,-1,0,1,2])
plt.ylabel("Sea Level (m)")

#lunar mixed
M2_Amplitude = 1.5
O1_Amplitude = 0.2
K1_Amplitude = 0.4
WaterLevels = M1_Amplitude*np.cos(np.radians(M1_Speed*TimesHours))
WaterLevels += M2_Amplitude*np.cos(np.radians(M2_Speed*TimesHours))
WaterLevels += O1_Amplitude*np.cos(np.radians(O1_Speed*TimesHours))
WaterLevels += K1_Amplitude*np.cos(np.radians(K1_Speed*TimesHours))
WaterLevels += N2_Amplitude*np.cos(np.radians(N2_Speed*TimesHours))
WaterLevels += K2_Amplitude*np.cos(np.radians(K2_Speed*TimesHours))
ax = plt.subplot(616)
plt.plot(TideTimes,WaterLevels,'k-')
ax.set_xticklabels([])
plt.yticks([-2,-1,0,1,2])
plt.ylabel("Sea Level (m)")
#WaterLevels += Mm_Amplitude*np.sin((2.*np.pi*TideTimes/(Mm_Period)))
#WaterLevels += S2_Amplitude*np.sin((2.*np.pi*TideTimes/(S2_Period)))

plt.xlabel("Time (hours)")
plt.show()