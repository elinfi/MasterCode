import kmedoids

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

from dissimilarity_matrix import DissimilarityMatrix


class KMedoids(DissimilarityMatrix):
    """Summary of class here.
    
    Longer class information...
    Longer class information...
    
    Attributes:
        list of everything with self.
        distmat: An array containing the parwise distances
    """
    
    def __init__(self, distmat):
        """Inits KMedoids class
        
        Args:
            distmat (ndarray):
                Distance matrix on square form
        """
        self.distmat = distmat
        
    def clusters(self, medoids, max_iter=100, init='random', random_state=None):
        """Calculates k-medoids clustering using fasterpam.
        
        Args:
            medoids (int):
                Number of clusters to find or existing medoids
            max_iter (int, optional):
                Maximum number of iterations
            init (str, optional):
                {'random', 'first', 'build'}, default = 'random'
                Initialization method
            randomstate (int, RandomState instance or None, optional):
                Random seed if no medoids are given
                
        Returns:
            clusters (KMedoidsResults):
                K-medoids clustering result
        """
        self.kmedoids_result = kmedoids.fasterpam(self.distmat, 
                                                  medoids, 
                                                  max_iter, 
                                                  init, 
                                                  random_state)
        return self.kmedoids_result
    
    def labels4plotting(self, n):
        """Convert cluster labels to 2d matrix for plotting.
        
        Args:
            n (int):
                Size of square matrix
        
        Return:
            labels_mat (ndarray):
                2D matrix with cluster labels
        """
        labels_list = self.kmedoids_result.labels # get cluster labels
        
        # initialize cluster matrix
        labels_mat = np.zeros((n, n))
        
        # fill matrix with cluster labels 
        upper_triag_idx = np.triu_indices(n)
        labels_mat[upper_triag_idx] = labels_list
        
        labels_mat = labels_mat + labels_mat.T
        labels_mat[np.diag_indices(n)] /= 2 # diagonal elements added twice
        
        return labels_mat
    
    def labels4plotting_nan(self, n, triu_nan_idx):
        """Convert cluster labels to 2d matrix for plotting. Insert removed NaN.
        
        Args:
            n (int):
                Size of square matrix
            triu_nan_idx (ndarray, dtype=bool):
                Boolean array with same length as the number of elements as
                the upper triangular matrix, with True for all NaN values and
                False elsewhere
        
        Return:
            labels_mat (ndarray):
                2D matrix with cluster labels
        """
        # index for upper triangular matrix
        upper_triag_idx = np.triu_indices(n)
        
        # get cluster labels
        labels_list = self.kmedoids_result.labels
        
        # empty list to fill with labels and original NaNs
        labels_nan_list = np.zeros_like(triu_nan_idx, dtype=float)

        labels_nan_list[~triu_nan_idx] = labels_list
        # labels_nan_list = np.where(~triu_nan_idx, labels_nan_list, np.nan)
        labels_nan_list[triu_nan_idx] = np.nan
        
        # initialize cluster matrix
        labels_mat = np.zeros((n, n))
        
        # fill matrix with cluster labels
        labels_mat[upper_triag_idx] = labels_nan_list
        
        labels_mat = labels_mat + labels_mat.T
        labels_mat[np.diag_indices(n)] /= 2 # diagonal elements were added twice
        
        self.labels_list = labels_list
        self.labels_nan_list = labels_nan_list
        
        return labels_mat