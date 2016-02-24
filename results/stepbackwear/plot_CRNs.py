# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:10:34 2016

@author: mhurst
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt

#open morphology file and read
FileName = "../driver_files/CRN_Model_Results.txt"
X, N = np.loadtxt(FileName,skiprows=1, unpack=True)
    
#setup figure
plt.figure(1,figsize=(6,3))
plt.plot(X,N,'k-')
plt.xlabel('Distance(m)')
plt.ylabel('$^{10}$Be concentration (atoms g$^{-1}$)')
plt.show()    