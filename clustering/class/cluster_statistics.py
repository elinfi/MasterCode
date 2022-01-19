"""Calculates statistics for the final clusters. 

Diagonal distance:
    Calculates the average distance from the diagonal for within each cluster.
    
Interactions mean:
    Calculates the mean number of interactions within each cluster.
    
Cluster mean:
    Calculates the average metric distance for each cluster.
"""

import numpy as np

class ClusterStatistics:
    def __init__(self, k, labels, X):
        """
        Args:
            k (int):
                Number of clusters.
            labels (ndarray):
                2D matrix with cluster labels.
            X (ndarray):
                n x 3 array containing the row indices in the first column,
                column indices in the second colum and the number of 
                interactions in the last column.
        """
        self.k = k
        self.labels = labels
        self.X = X
        
    def diagonal_dist(self, data):
        """Calculates the average distance from the diagonal.
        
        Args:
            data (ndarray):
                n x 3 array containing the row indices in the first column,
                column indices in the second colum and the number of 
                interactions in the last column for one cluster.
        Returns:
            diag_mean (float):
                The average distance from the diagonal for one cluster.
        """
        dist = abs(data[:, 0] - data[:, 1])
        diag_mean = np.mean(dist)
        return diag_mean
    
    def interactions_dist(self, data):
        """Calculates the mean of interactions.
        
        Args:
            data (ndarray):
                n x 3 array containing the row indices in the first column,
                column indices in the second colum and the number of 
                interactions in the last column for one cluster.
        Returns:
            interactions_mean (float):
                The average number of interactions for one cluster.
        """
        interactions_mean = np.mean(data[:, 2])
        return interactions_mean
        
    def cluster_mean(self, metric):
        """Calculates the average metric distance for each clusters.
        
        Args:
            metric (str) {'row_col_dist'}:
                Which metric to use when calculating the dissimilarity matrix.
        Returns:
            mean_list (ndarray):
                Array with the average distance for each cluster.
        """
        metric_func = getattr(self, metric)
        
        mean_list = np.zeros(self.k)
        
        for i in range(self.k):
            data = self.X[self.labels == i]
            mean_list[i] = metric_func(data)
            
        return mean_list