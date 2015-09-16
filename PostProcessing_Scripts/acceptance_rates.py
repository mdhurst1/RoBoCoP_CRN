# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 10:21:28 2015

Determine the acceptance rate for each scenario

@author: mhurst
"""

#IMPORT MODULES
import numpy as np
import matplotlib.pyplot as plt

Filenames = ["HG_single_chain2.dout","BH_single_chain2.dout",
             "HG_linear_chain2.dout","BH_linear_chain2.dout",
             "HG_change_chain2.dout","BH_change_chain2.dout"]

for Filename in Filenames:
    
    #Load the chain file
    i, RetreatRate1, BeachWidth, NewLike, LastLike, NAccepted, NRejected  = np.loadtxt(Filename, unpack=True, skiprows=2, usecols=(0,1,4,6,7,8,9))

    # Get accepted parameters for plotting convergence
    N = 0
    for j in range(0,len(i)):
        if NAccepted[j] > N:
            N += 1
            

    print Filename + "\n\t" + \
            str(N) + "/" + str(len(i)) + " accepted" + "\n\t" + \
            str(100.*float(N)/float(len(i))) + "% accepted"
    