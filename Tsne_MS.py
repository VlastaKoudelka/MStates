# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 13:32:48 2016

@author: vlastimilo
"""

import tsne
import matplotlib.pyplot as plt
import numpy as np

perplexity = 5
initdims = 19
outdims = 2

file_path = '/home/vlastimilo/NUDZ_Data/Microstates/'

#file_path = 'D:\MS_data\\'
file_name = 'BAP-Brisova-avg-19ele-GFPpeaks.txt'

GFP_data = np.loadtxt(file_path+file_name)
#GFP_data.astype(int)
GFP_data = GFP_data/np.amax(GFP_data)
#GFP_data = GFP_data[450:500]

[mapped,C] = tsne.tsne(GFP_data, outdims, initdims, perplexity)

plt.figure
plt.scatter(mapped[:,0],mapped[:,1])
plt.title('t-SNE GFP_peaks 5 perplexity')

f_name = 't-SNE GFP_peaks_perp_5.jpeg'
plt.savefig(f_name, dpi=300, facecolor='w', edgecolor='w',
orientation='portrait', papertype=None, format=None,
transparent=False, bbox_inches=None, pad_inches=0.1,
frameon=None)    