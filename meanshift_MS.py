# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 10:26:40 2016
Mean-shift
@author: vlastimilo
"""

import matplotlib.pyplot as plt
import numpy as np
import sklearn.cluster as cluster
from scipy.spatial.distance import pdist as pdist

file_path = '/home/vlastimilo/NUDZ_Data/Microstates/' #LINUX
#file_path = 'D:\MS_data\\'  #WINDOWS
file_name = 'BAP-Brisova-avg-19ele-GFPpeaks.txt'

GFP_data = np.loadtxt(file_path+file_name)

GFP_dist_mat = pdist(GFP_data)
GFP_mean_dist = np.average(GFP_dist_mat)
GFP_mean_dist = GFP_mean_dist

GFP_obj = cluster.MeanShift(bandwidth=GFP_mean_dist, seeds=None, bin_seeding=False, min_bin_freq=1, cluster_all=True, n_jobs=4)
GFP_labels = GFP_obj.fit_predict(GFP_data)
print "A number of clusters:" + str(np.max(GFP_labels))

kernels = np.divide(GFP_mean_dist,np.linspace(1,6,20))

GFP_n_clust = kernels
for i,width in enumerate(kernels):
    GFP_obj = cluster.MeanShift(bandwidth=width, seeds=None, bin_seeding=False, min_bin_freq=1, cluster_all=True, n_jobs=4)
    GFP_labels = GFP_obj.fit_predict(GFP_data)
    GFP_n_clust[i] = np.max(GFP_labels)
    
plt.figure
plt.plot(GFP_n_clust)