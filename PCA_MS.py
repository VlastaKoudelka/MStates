# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 13:48:56 2016

@author: vlastimilo
"""

import matplotlib.pyplot as plt
import numpy as np
import sklearn.cluster as cluster
from sklearn.decomposition import PCA 
from scipy.spatial.distance import pdist as pdist

file_path = '/home/vlastimilo/NUDZ_Data/Microstates/' #LINUX
#file_path = 'D:\MS_data\\'  #WINDOWS
file_name = 'BAP-Brisova-avg-19ele-GFPpeaks.txt'

GFP_data = np.loadtxt(file_path+file_name)
pca_obj = PCA()
pca_obj.fit(GFP_data)
print(pca_obj.explained_variance_ratio_)
plt.figure
plt.plot(pca_obj.explained_variance_ratio_)
plt.xlabel('number of components')
plt.ylabel('explained variance')