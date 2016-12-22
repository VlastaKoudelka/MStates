# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 10:33:02 2016

@author: vlastimilo
"""

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

neighborhood = np.divide(GFP_mean_dist,np.linspace(1,6,20))


GFP_n_clust = np.zeros((7,len(neighborhood)))
plt.figure
for j,min_sapl in enumerate(np.arange(7)+3):  
    for i,width in enumerate(neighborhood):
        print i
        
        dbscan_obj = cluster.DBSCAN(eps=width, min_samples=min_sapl)    
        dbscan_labels = dbscan_obj.fit_predict(GFP_data)
        
        GFP_n_clust[j,i] = np.max(dbscan_labels)
    
    plt.plot(np.linspace(1,6,20),GFP_n_clust[j],label = 'no_of_neighbors: '+ str(min_sapl))



plt.xlabel('mean distance division (distance in one neighborhood)')
plt.ylabel('no. of clusters estimate')
plt.title('DBSCAN algorithm cluster estimation')
plt.legend()

f_name = 'dbscan_no_clst.jpeg'
plt.savefig(f_name, dpi=300, facecolor='w', edgecolor='w',
orientation='portrait', papertype=None, format=None,
transparent=False, bbox_inches=None, pad_inches=0.1,
frameon=None)  