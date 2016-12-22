# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 13:31:07 2016
k-means
@author: vlastimilo
"""

import matplotlib.pyplot as plt
import numpy as np
import sklearn.cluster as cluster
from sklearn import metrics

file_path = '/home/vlastimilo/NUDZ_Data/Microstates/' #LINUX
#file_path = 'D:\MS_data\\'  #WINDOWS
file_name = 'BAP-Brisova-avg-19ele-GFPpeaks.txt'

GFP_data = np.loadtxt(file_path+file_name)

crit_kmeans = np.zeros((3,9))
for i,n_clst in enumerate(np.arange(2,11,1)):
    obj_kmeans = cluster.KMeans(n_clst)
    labels_kmeans = obj_kmeans.fit_predict(GFP_data)
    print i
    obj_aglo = cluster.AgglomerativeClustering(n_clst)
    labels_aglo = obj_aglo.fit_predict(GFP_data)
    
#    obj_spect = cluster.SpectralClustering(n_clst)
#    labels_spect = obj_spect.fit_predict(GFP_data)
    
    crit_kmeans[0,i] = metrics.silhouette_score(GFP_data,labels_kmeans)
    crit_kmeans[1,i] = metrics.silhouette_score(GFP_data,labels_aglo)
#    crit_kmeans[2,i] = metrics.silhouette_score(GFP_data,labels_aglo)

plt.figure
x =np.arange(2,11,1)
plt.plot(x,crit_kmeans[0],label='Kmeans')
plt.plot(x,crit_kmeans[1],label='Agglomerative')
#plt.plot(x,crit_kmeans[2],label='Spectral')
plt.title('Informed algorithm clustering')
plt.xlabel('number of clusters')
plt.ylabel('silhouette criterion')
plt.legend()
plt.show    

f_name = 'Informed_algorithms.jpeg'
plt.savefig(f_name, dpi=300, facecolor='w', edgecolor='w',
orientation='portrait', papertype=None, format=None,
transparent=False, bbox_inches=None, pad_inches=0.1,
frameon=None)  