# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 13:32:48 2016

@author: vlastimilo
"""

import tsne
import matplotlib.pyplot as plt
import numpy as np

perplexity = 20
initdims = 19
outdims = 2

file_path = '/home/vlastimilo/NUDZ_Data/Microstates/'
file_name = 'BAP-Brisova-avg-19ele-GFPpeaks.txt'

GFP_data = np.loadtxt(file_path+file_name)

[mapped,C] = tsne.tsne(GFP_data, outdims, initdims, perplexity)

plt.figure
plt.scatter(mapped[:,0],mapped[:,1])