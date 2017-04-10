# -*- coding: utf-8 -*-
"""
Created on Wed Apr 05 11:40:51 2017

Milnkovitch cyclicity of sea level

@author: Martin Hurst
"""

# import modules
import numpy as np
import matplotlib.pyplot as plt

# setup control parmeters
Time = np.arange(-700000,10000)

#Obliquity (tilt)
#22.1 - 24.5 degrees variation
Offset = np.random.rand()
O_Mean = 0
O_Amplitude = 5
O_WaveLength = 41000
O_Signal = O_Mean+O_Amplitude*np.cos(2.*np.pi*(Time+Offset*O_WaveLength)/O_WaveLength)
#plt.plot(Time,WaterLevel)
#plt.show()

#Precision (wobble)
Offset = np.random.rand()
P_Mean = 0
P_Amplitude = 8
P_WaveLength = 26000
P_Signal = P_Mean+P_Amplitude*np.cos(2.*np.pi*(Time+Offset*P_WaveLength)/P_WaveLength)

#Eccentricity
# 413,000 years, and changes the eccentricity by ±0.012. 
# Two other components have periods of 95,000 and 125,000
E_Mean = 0
E_Amplitude = 90
E_WaveLength = 120000
Offset = 10000
E_Signal = (E_Amplitude/np.pi)*np.arctan(1./np.tan(np.pi*(Time+Offset)/E_WaveLength))
E_Signal += E_Mean+0.1*E_Amplitude*np.sin(2.*np.pi*Time/E_WaveLength)
#E_Signal *= 1


#Total signal
Total_Signal = O_Signal+P_Signal+E_Signal

print np.max(Total_Signal)-np.min(Total_Signal)

#make the plot
plt.figure(1,figsize=(6.6,6.6))
Time /= 1000
plt.subplot(411)
plt.plot(Time,O_Signal,'g',lw=1.5)
plt.xticks([])
plt.ylim(np.min(O_Signal)+0.5*np.min(O_Signal),np.max(O_Signal)+0.5*np.max(O_Signal))
plt.ylabel("contribution (m)")
plt.text(-700,np.max(O_Signal)+0.1*np.max(O_Signal),'Obliquity')
plt.subplot(412)
plt.plot(Time,P_Signal,'b',lw=1.5)
plt.xticks([])
plt.ylim(np.min(P_Signal)+0.5*np.min(P_Signal),np.max(P_Signal)+0.5*np.max(P_Signal))
plt.ylabel("contribution (m)")
plt.text(-700,np.max(P_Signal)+0.1*np.max(P_Signal),'Precession')
plt.subplot(413)
plt.plot(Time,E_Signal,'r',lw=1.5)
plt.xticks([])
plt.ylim(np.min(E_Signal)+0.5*np.min(E_Signal),np.max(E_Signal)+0.5*np.max(E_Signal))
plt.ylabel("contribution (m)")
plt.text(-700,np.max(E_Signal)+0.1*np.max(E_Signal),'Eccentricity')
plt.subplot(414)
plt.plot(Time,Total_Signal,'k-',lw=3)
plt.ylim(np.min(Total_Signal)+0.5*np.min(Total_Signal),np.max(Total_Signal)+0.5*np.max(Total_Signal))
plt.xlabel('Time ky')
plt.ylabel('Sea level (m)')
plt.text(-700,np.max(Total_Signal)+0.1*np.max(Total_Signal),'Total Signal')
plt.tight_layout()
plt.show()