"""
####################
clustering.py
####################
Code supporting efforts to optimize and automate external clustering libraries for gFunc purposes.
"""

import Pycluster as pClstr
import mlpy

import numpy as np
import matplotlib.pyplot as plt

def __generate_random_data_centers(size1=200,size2=300,size3=200):
    """
    *DOES:*
        * generates 3 sets of data containing randomized members centered around known coords (x1,x2,x3)
    *RETURNS:*
        * tuple(data,x1,x2,x3)
    """
    
    mean1, cov1, n1 = [1, 5], [[1,1],[1,2]], size1 # 200 points, mean=(1,5)
    x1 = np.random.multivariate_normal(mean1, cov1, n1)
    mean2, cov2, n2 = [2.5, 2.5], [[1,0],[0,1]], size2 # 300 points, mean=(2.5,2.5)
    x2 = np.random.multivariate_normal(mean2, cov2, n2)
    mean3, cov3, n3 = [5, 8], [[0.5,0],[0,0.5]], size3 # 200 points, mean=(5,8)
    x3 = np.random.multivariate_normal(mean3, cov3, n3)
    data = np.concatenate((x1, x2, x3), axis=0) # concatenate the samples
    membership = np.array([0]*size1 + [1]*size2 + [2]*size3)
    return (data,x1,x2,x3,membership)

def plot_centers_and_points(data,clusters,means,truth=None):
    """
    """
    fig = plt.figure(1)
    plot1 = plt.scatter(data[:,0], data[:,1], c=clusters, s=80, alpha=.5, marker='s')
    if truth is not None:
        plotT = plt.scatter(data[:,0], data[:,1], c=truth, s=30, alpha=.5, marker='o') # plot KNOWN memberships
    plot2 = plt.scatter(means[:,0], means[:,1], c=np.unique(clusters), s=128, marker='d') # plot the means
    plt.show()
    