# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 13:32:48 2016

@author: vlastimilo
"""

import tsne
import matplotlib.pyplot as plt
import numpy as np

perplexity = 50
initdims = 19
outdims = 2

#file_path = '/home/vlastimilo/NUDZ_Data/Microstates/'

file_path = 'D:\MS_data\\'
file_name = 'BAP-Brisova-avg-19ele-GFPpeaks.txt'

GFP_data = np.loadtxt(file_path+file_name)
#GFP_data.astype(int)
GFP_data = GFP_data/np.amax(GFP_data)
#GFP_data = GFP_data[450:500]

[mapped,C] = tsne.tsne(GFP_data, outdims, initdims, perplexity)

plt.figure
plt.scatter(mapped[:,0],mapped[:,1])