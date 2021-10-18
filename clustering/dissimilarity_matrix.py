import kmedoids

import numpy as np
import numpy.ma as ma

from numba import jit
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import MinMaxScaler



class DissimilarityMatrix:
    def __init__(self, data):
        """Creates a dissimilarity matrix.

        Creates a dissimilarity based on a Hi-C data matrix with a given
        distance metric.

        Args:
            data (ndarray):
                Square symmetric matrix containing the Hi-C data to cluster
            metric {'euclidean', 'test_euclidean'}:
                Which metric to use when calculating the dissimilarity matrix
        Attributes:
            n (int):
                Length of data matrix
        """
        self.data = data
        self.n = data.shape[0]
        
    def interactions_dist(self, u, v):
        """Calculates the difference in number of interactions.
        
        Calculates the difference in number of interactions given by
            dist = |u_z - v_z|
        where z is the number of interactions.
        
        Args:
            u (ndarray):
                Array containing one data point.
            v (ndarray):
                Array containing one data point.
        Result:
            dist (float):
                Difference in number of interactions between u and v.
        """
        dist = abs(u[2] - v[2])
        return dist
        
    def diagonal_dist(self, u, v):
        """Calculates the scaled difference in distance from the main diagonal.
        
        Calculates the scaled  difference in distance from the main diagonal 
        between u and v given by
            dist = | |u_i - v_i| - |u_j - v_j| | / n,
        where i is the row index, j is the column index and n is the maximum 
        distance from the main diagonal. The maximum diagonal distance equals
        the length of the matrix when square.

        Args:
            u (ndarray):
                Array containing one data point.
            v (ndarray):
                Array containing one data point.
        Result:
            dist (float):
                Difference in distance from the main diagonal between u and v.
        """
        d1 = abs(u[0] - u[1])
        d2 = abs(v[0] - v[1])
        
        dist = abs(d1 - d2)/self.n
        return dist
    
    def manhattan(self, u, v):
        """Calculates the scaled manahattan distance between u and v.
        
        Calculates the scaled manhattan distance between u and v given by
            dist = (|u_i - v_i| + |u_j - v_j|)/2n,
        where i is the row index,j is the column index and 2n is the maximum 
        manhattan distance.
        
        Args:
            u (ndarray):
                Array containing one data point.
            v (ndarray):
                Array containing one data point.
        Result:
            dist (float):
                Manahattan distance between u and v.
        """
        dist = np.sum(u - v)/(2*self.n)
        return dist

    def transform_data_nan(self, scaler=None):
        """Extract the upper triangular of the data matrix.

        Return:
            X (ndarray):
                n x 3 array containing the bin1 and bin2 in the first two
                columns and the number of interactions in the last column
        """
        #assert log_scale in [None, 'log2', 'log10', 'MinMax'], 'Only None, 2 or 10 are valid log scales'
            
        n = self.n
        self.triu_idx = np.triu_indices(n)
        self.nan_idx = np.isnan(self.data)
        self.triu_nan_idx = self.nan_idx[self.triu_idx]

        # create 2d coordinate system in 1d of the upper triangular matrix
        xy = np.ones((n, n))*np.linspace(0, n-1, n) 
        x = xy.T[self.triu_idx]
        x = x[~self.triu_nan_idx]
        y = xy[self.triu_idx]
        y = y[~self.triu_nan_idx]
        
        # extract the upper triangular matrix of data into a flat array      
        z = self.data[self.triu_idx]
        z = z[~self.triu_nan_idx]
        if scaler:
            scaler.fit(z.reshape(1, -1))
            z_scaled = scaler.transform(z.reshape(1, -1))
            z = z_scaled[0, :]
        #if log_scale == 2:
        #    z = np.log2(z)
        #elif log_scale == 10:
        #    z = np.log10(z)
            

        # create matrix with x, y, z as columns
        X = np.array([x, y, z]).T

        return X
 
    def transform_data(self, scaler=None):
        """Extract the upper triangular of the data matrix.

        Return:
            X (ndarray):
                n x 3 array containing the bin1 and bin2 in the first two
                columns and the number of interactions in the last column
        """
        #assert log_scale in [None, 'log2', 'log10', 'MinMax'], 'Only None, 2 or 10 are valid log scales'
            
        n = self.n
        triu_idx = np.triu_indices(n)

        # create 2d coordinate system in 1d of the upper triangular matrix
        xy = np.ones((n, n))*np.linspace(0, n-1, n)
        x = xy.T[triu_idx]
        y = xy[triu_idx]
        
        # extract the upper triangular matrix of data into a flat array      
        z = self.data[triu_idx]
        if scaler:
            scaler.fit(z.reshape(1, -1))
            z_scaled = scaler.transform(z.reshape(1, -1))
            z = z_scaled[0, :]
        #if log_scale == 2:
        #    z = np.log2(z)
        #elif log_scale == 10:
        #    z = np.log10(z)
            

        # create matrix with x, y, z as columns
        X = np.array([x, y, z]).T

        return X
          
    def scipy_dist(self, metric, scaler=None, remove_nan=True, **kwargs):
        """
        Computes the distance matrix using scipy.pdist with a given metric.
        
        Args:
            metric {'row_col_dist'}:
                Which metric to use when calculating the dissimilarity matrix
        Return:
            distmat (array):
                Distance matrix on square form
        """
        # transform data
        if remove_nan:
            X = self.transform_data_nan(scaler)
        else:
            X = self.transform_data(scaler)
        
        if hasattr(self, metric):
            # call metric function
            metric = getattr(self, metric)
            
        # calculate the pairwise distances between data point and convert it
        # to square distance matrix
        distmat = squareform(pdist(X, metric=metric, **kwargs))

        return distmat
    
    def sklearn_dist(self, metric, scaler=None, **kwargs):
        """
        Computes the distance matrix using sklearn.pairwise_distances.
        
        Args:
            metric {'euclidean', 'test_euclidean'}:
                Which metric to use when calculating the dissimilarity matrix
        Returns:
            distmat (array):
                Distance matrix on square form
        """
        # transform data
        X = self.transform_data(scaler)        
        
        if hasattr(self, metric):
            # call metric function
            metric = getattr(self, metric)

        # calculate the pairwise distances between data point and convert it
        # to square distance matrix
        distmat = pairwise_distances(X, metric=metric, force_all_finite=False, 
                                     n_jobs=-1, **kwargs)
        
        return distmat


if __name__ == '__main__':
    test = np.array([[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4]])
    distmat = DissimilarityMatrix(test)
    print(distmat.scipy_dist('minkowski_', log_scale=None, k=2))
