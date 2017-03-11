# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 02:53:09 2017

@author: martin
"""

import numpy as np
import matplotlib.pyplot as plt

XMax = 10
x = np.arange(0.01,XMax,0.01)
Theta = 0
m = 0.0805*XMax+0.25
sigma = 0.5
print np.exp(m-sigma**2.)

P = (np.exp(-((np.log(x-Theta)-m)**2.)/(2*sigma**2.))) / ((x-Theta)*sigma*np.sqrt(2*np.pi))
P = P/np.max(P)
plt.plot(x,P)
plt.show()